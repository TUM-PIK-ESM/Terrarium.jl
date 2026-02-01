"""
    $TYPEDEF

Represents an explicit formulation of the two-phase heat conduction operator in 1D:

```math
\\frac{\\partial U(T,\\phi)}{\\partial t} = \\nabla \\cdot \\kappa(T)\\nabla_x T(x,t)
```
where \$T\$ is temperature [K], \$U\$ is internal energy [J m⁻³], and \$\\kappa\$ is the thermal conductivity [W m K⁻¹].
"""
@kwdef struct ExplicitTwoPhaseHeatConduction <: AbstractHeatOperator end

"""
    $TYPEDEF

Standard implementation of the soil energy balance accounting for freezing and thawing of pore water/ice.
The `closure` field represents the temperature-energy closure \$U(T,\\phi)\$ which relates temperature to internal
energy via an arbitrary set of additional parameters \$\\phi\$ which are determined by the model configuration.

Properties:
$TYPEDFIELDS
"""
struct SoilEnergyBalance{
    NF,
    HeatOperator<:AbstractHeatOperator,
    EnergyClosure<:AbstractSoilEnergyClosure,
    ThermalProps<:SoilThermalProperties{NF}
} <: AbstractSoilEnergyBalance{NF}
    "Heat transport operator"
    operator::HeatOperator

    "Closure relating energy and temperature"
    closure::EnergyClosure

    "Soil thermal properties"
    thermal_properties::ThermalProps
end

SoilEnergyBalance(
    ::Type{NF};
    operator::AbstractHeatOperator = ExplicitTwoPhaseHeatConduction(),
    closure::AbstractSoilEnergyClosure = SoilEnergyTemperatureClosure(),
    thermal_properties::SoilThermalProperties{NF} = SoilThermalProperties(NF),
) where {NF} = SoilEnergyBalance(operator, closure, thermal_properties)

variables(energy::SoilEnergyBalance) = (
    prognostic(:internal_energy, XYZ(); closure=energy.closure, units=u"J/m^3", desc="Internal energy of the soil volume, including both latent and sensible components"),
    auxiliary(:ground_temperature, XY(), ground_temperature, energy, units=u"°C", desc="Temperature of the uppermost ground or soil grid cell in °C"),
)

# Field constructor for ground_temperature that returns a view of the uppermost soil layer
function ground_temperature(energy::SoilEnergyBalance, grid, clock, fields)
    fgrid = get_field_grid(grid)
    # Use uppermost soil layer as ground temperature
    # TODO: Revisit this if/when we extend the vertical layers to include snow and canopy
    return @view fields.temperature[:, :, fgrid.Nz]
end

get_closure(energy::SoilEnergyBalance) = energy.closure

function initialize!(
    state, grid,
    energy::SoilEnergyBalance,
    ground::AbstractGround,
    constants::PhysicalConstants,
    args...
)
    # Initialize by evaluating inverse closure (temperature -> energy)
    # Note that this assumes the temperature state to have already been initialized!
    invclosure!(state, grid, energy, ground, constants)
    return nothing
end

compute_auxiliary!(state, grid, energy::SoilEnergyBalance, args...) = nothing

function compute_tendencies!(
    state, grid,
    energy::SoilEnergyBalance,
    ground::AbstractGround,
    args...
)
    kernel_args = (
        energy,
        get_hydrology(ground),
        get_stratigraphy(ground),
        get_biogeochemistry(ground)
    )
    fields = get_fields(state, kernel_args...)
    tendencies = tendency_fields(state, energy)
    launch!(grid, XYZ, compute_energy_tendency_kernel!, tendencies, fields, kernel_args)
    return nothing
end

# Kernel functions

@inline function compute_energy_tendencies!(
    tendencies, i, j, k, grid, fields,
    energy::SoilEnergyBalance,
    args...
)
    tendencies.internal_energy[i, j, k] += compute_energy_tendency(i, j, k, grid, fields, energy, args...)
end

@inline function compute_energy_tendency(
    i, j, k, grid, fields,
    energy::SoilEnergyBalance{NF, <:ExplicitTwoPhaseHeatConduction},
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    # Operators require the underlying Oceananigans grid
    field_grid = get_field_grid(grid)

    # Divergence of heat fluxes
    ∂U∂t = -∂zᵃᵃᶜ(i, j, k, field_grid, diffusive_heat_flux, fields, energy, hydrology, strat, bgc)
    return ∂U∂t
end

@inline function compute_thermal_conductivity(
    i, j, k, grid, fields,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    soil = soil_volume(i, j, k, grid, fields, strat, hydrology, bgc)
    return compute_thermal_conductivity(energy.thermal_properties, soil)
end

# Diffusive heat flux term passed to ∂z operator
@inline function diffusive_heat_flux(i, j, k, grid, fields, args...)
    # Get temperature field
    T = fields.temperature
    # Compute and ∂U∂tinterpolate conductivity to grid cell faces
    κ = ℑzᵃᵃᶠ(i, j, k, grid, compute_thermal_conductivity, fields, args...)
    # Fourier's law: q = -κ ∂T/∂z
    q = -κ * ∂zᵃᵃᶠ(i, j, k, grid, T)
    return q
end

# Kernels

@kernel function compute_tendencies_kernel!(tendencies, grid, fields, energy::SoilEnergyBalance, args)
    i, j, k = @index(Global, NTuple)
    compute_energy_tendencies!(tendencies, i, j, k, grid, fields, energy, args...)
end
