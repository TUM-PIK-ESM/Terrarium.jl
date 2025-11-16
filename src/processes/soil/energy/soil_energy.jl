abstract type AbstractHeatOperator <: AbstractOperator end

"""
    $TYPEDEF

Represents an explicit formulation of the two-phase heat conduction operator in 1D:

```math
\\frac{\\partial U(T,\\phi)}{\\partial t} = \\nabla \\cdot \\kappa(T)\\nabla_x T(x,t)
```
where \$T\$ is temperature [K], \$U\$ is internal energy [J m⁻³], and \$\\kappa\$ is the thermal conductivity [W m K⁻¹].
The `closure` field represents the temperature-energy closure \$U(T,\\phi)\$ which relates temperature to internal
energy via an arbitrary set of additional parameters \$\\phi\$ which are determined by the model configuration.
"""
@kwdef struct ExplicitTwoPhaseHeatConduction{ET} <: AbstractHeatOperator
    closure::ET = TemperatureEnergyClosure()
end

get_closure(op::ExplicitTwoPhaseHeatConduction) = op.closure

"""
    $TYPEDEF

Standard implementation of the soil energy balance accounting for freezing and thawing of pore water/ice.

Properties:
$TYPEDFIELDS
"""
struct SoilEnergyBalance{
    NF,
    HeatOperator<:AbstractHeatOperator,
    FC<:FreezeCurve,
    ThermalProps<:SoilThermalProperties{NF}
} <: AbstractSoilEnergyBalance{NF}
    "Heat transport operator"
    operator::HeatOperator

    "Soil thermal properties"
    thermal_properties::ThermalProps

    "Freeze curve type constructor"
    freezecurve::Val{FC}
end

SoilEnergyBalance(
    ::Type{NF};
    operator::AbstractHeatOperator = ExplicitTwoPhaseHeatConduction(),
    thermal_properties::SoilThermalProperties{NF} = SoilThermalProperties(NF),
    fctype::Type{<:FreezeCurve} = FreeWater
) where {NF} = SoilEnergyBalance(operator, thermal_properties, Val{fctype}())

freezecurve(
    ::SoilEnergyBalance{NF, OP, FreeWater},
    ::AbstractSoilHydrology
) where {NF, OP} = FreeWater()

variables(energy::SoilEnergyBalance) = (
    prognostic(:temperature, XYZ(), closure=get_closure(energy.operator), units=u"°C", desc="Temperature of the grid cell in °C"),
    auxiliary(:saturation_water_ice, XYZ(), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
    auxiliary(:liquid_water_fraction, XYZ(), domain=UnitInterval(), desc="Fraction of unfrozen water in the pore space"),
)


# evaluate closure (temperature -> energy) on initialize!
function initialize!(state, model, energy::SoilEnergyBalance)
    closure!(state, model, state.closures.temperature)
    return nothing
end

compute_auxiliary!(state, model, energy::SoilEnergyBalance) = nothing

function compute_tendencies!(state, model, energy::SoilEnergyBalance)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    launch!(grid, :xyz, compute_energy_tendency!, state, grid, energy, hydrology, strat, bgc)
    return nothing
end

@kernel function compute_energy_tendency!(
    state,
    grid,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
)
    idx = @index(Global, NTuple)
    state.tendencies.internal_energy[idx...] += energy_tendency(idx, state, grid, energy, hydrology, strat, bgc)
end

@inline function energy_tendency(
    idx, state, grid,
    energy::SoilEnergyBalance{NF, <:ExplicitTwoPhaseHeatConduction},
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
) where {NF}
    i, j, k = idx
    # operators require the underlying Oceananigans grid
    field_grid = get_field_grid(grid)

    # Divergence of heat fluxes
    dUdt = -∂zᵃᵃᶜ(i, j, k, field_grid, diffusive_heat_flux, state, energy, hydrology, strat, bgc)
    return dUdt
end

@inline function thermalconductivity(
    i, j, k, grid, state,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
)
    idx = (i, j, k)
    soil = soil_composition(idx, state, strat, hydrology, bgc)
    return thermalconductivity(energy.thermal_properties, soil)
end

# Diffusive heat flux term passed to ∂z operator
@inline function diffusive_heat_flux(i, j, k, grid, state, args...)
    # Get temperature field
    T = state.temperature
    # Compute and interpolate conductivity to grid cell faces
    κ = ℑzᵃᵃᶠ(i, j, k, grid, thermalconductivity, state, args...)
    # Fourier's law: q = -κ ∂T/∂z
    q = -κ * ∂zᵃᵃᶠ(i, j, k, grid, T)
    return q
end
