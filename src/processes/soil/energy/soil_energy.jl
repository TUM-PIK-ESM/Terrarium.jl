abstract type AbstractHeatOperator end

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
    closure::ET = EnergyTemperatureClosure()
end

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

get_closure(energy::SoilEnergyBalance) = energy.operator.closure

variables(energy::SoilEnergyBalance) = (
    prognostic(:internal_energy, XYZ(); closure=get_closure(energy), units=u"J/m^3", desc="Internal energy of the soil volume, including both latent and sensible components"),
    auxiliary(:liquid_water_fraction, XYZ(), domain=UnitInterval(), desc="Fraction of unfrozen water in the pore space"),
)

# evaluate inverse closure (temperature -> energy) on initialize!
function initialize!(state, model, energy::SoilEnergyBalance)
    invclosure!(state, model, get_closure(energy))
    return nothing
end

compute_auxiliary!(state, model, energy::SoilEnergyBalance) = nothing

function compute_tendencies!(state, model, energy::SoilEnergyBalance)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)

    model_compents = getter(model)
    launch!(grid, :xyz, compute_energy_tendency!, state, model_compents...)
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
    i, j, k = @index(Global, NTuple)
    state.tendencies.internal_energy[i, j, k] += energy_tendency(i, j, k, state, grid, energy, hydrology, strat, bgc)
end

@inline function energy_tendency(
    i, j, k, state, grid,
    energy::SoilEnergyBalance{NF, <:ExplicitTwoPhaseHeatConduction},
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
) where {NF}
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
    soil = soil_composition(i, j, k, state, strat, hydrology, bgc)
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

"""
    EnergyTemperatureClosure

Defines the constitutive relationship between temperature and the internal energy, U, of the system, i.e:

```math
U(T) = T\\times C(T) - L_f \\theta_{wi} (1 - F(T))
```
where T is temperature, C(T) is the temperature-dependent heat capacity, L_f is the
volumetric latent heat of fusion, and F(T) is the constitutive relation between temperature
and the unfrozen fraction of pore water. Note that, under this formulation, zero energy corresponds to
0°C with no ice, i.e. all pore water fully thawed.

The closure relation is defined as being a mapping from the conserved quantity (energy) to the continuous
quantity (temperature), i.e. the inverse of U(T).
"""
struct EnergyTemperatureClosure <: AbstractClosureRelation end

closurevar(::EnergyTemperatureClosure) = auxiliary(:temperature, XYZ(), units=u"°C", desc="Temperature of the soil volume in °C")

function closure!(state, model::AbstractSoilModel, ::EnergyTemperatureClosure)
    grid = get_grid(model)
    energy = get_soil_energy_balance(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    constants = get_constants(model)
    launch!(grid, :xyz, energy_to_temperature!, state, energy, hydrology, strat, bgc, constants)
    return nothing
end

function invclosure!(state, model::AbstractSoilModel, ::EnergyTemperatureClosure)
    grid = get_grid(model)
    energy = get_soil_energy_balance(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    constants = get_constants(model)
    launch!(grid, :xyz, temperature_to_energy!, state, energy, hydrology, strat, bgc, constants)
    return nothing
end

@kernel function temperature_to_energy!(
    state,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
)
    i, j, k = @index(Global, NTuple)
    fc = freezecurve(energy, hydrology)
    temperature_to_energy!(i, j, k, state, fc, energy, hydrology, strat, bgc, constants)
end

@inline function temperature_to_energy!(
    i, j, k, state, ::FreeWater,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
)
    T = state.temperature[i, j, k] # assumed given
    L = constants.ρw*constants.Lsl
    por = porosity(i, j, k, state, hydrology, strat, bgc)
    sat = state.saturation_water_ice[i, j, k]
    # calculate unfrozen water content from temperature
    # N.B. For the free water freeze curve, the mapping from temperature to unfrozen water content
    # within the phase change region is indeterminate since it is assumed that T = 0. As such, we
    # have to assume here that the liquid water fraction is zero if T < 0 and one otherwise. This method
    # should therefore only be used for initialization and should **not** be involved in the calculation
    # of tendencies.
    liq = state.liquid_water_fraction[i, j, k] = ifelse(
        T >= zero(T),
        one(sat),
        zero(sat),
    )
    soil = soil_composition(i, j, k, state, strat, hydrology, bgc)
    C = heatcapacity(energy.thermal_properties, soil)
    # compute energy from temperature, heat capacity, and ice fraction
    U = state.internal_energy[i, j, k] = T*C - L*sat*por*(1 - liq)
    return U
end

@kernel function energy_to_temperature!(
    state,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
)
    i, j, k = @index(Global, NTuple)
    fc = freezecurve(energy, hydrology)
    energy_to_temperature!(i, j, k, state, fc, energy, hydrology, strat, bgc, constants)
end

@inline function energy_to_temperature!(
    i, j, k, state, fc::FreeWater,
    energy::SoilEnergyBalance{NF},
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
) where {NF}

    U = state.internal_energy[i, j, k] # assumed given
    L = constants.ρw*constants.Lsl
    por = porosity(i, j, k, state, hydrology, strat, bgc)
    sat = state.saturation_water_ice[i, j, k]
    Lθ = L*sat*por
    # calculate unfrozen water content
    state.liquid_water_fraction[i, j, k] = liquid_water_fraction(fc, U, Lθ, sat)
    soil = soil_composition(i, j, k, state, strat, hydrology, bgc)
    C = heatcapacity(energy.thermal_properties, soil)
    # calculate temperature from internal energy and liquid water fraction
    T = state.temperature[i, j, k] = energy_to_temperature(fc, U, Lθ, C)
    return T
end

"""
Calculate the unfrozen water content from the given internal energy, latent heat content, and saturation.
"""
@inline function liquid_water_fraction(::FreeWater, U::NF, Lθ::NF, sat::NF) where {NF}
    return if U >= zero(U)
        # Case 1: U ≥ Lθ -> thawed
        one(sat)
    else
        # Case 2a: -Lθ ≤ U ≤ 0 -> phase change
        # Case 2b: U < -Lθ -> frozen (zero)
        (U >= -Lθ)*(one(sat) - safediv(U, -Lθ))
    end
end

"""
Calculate the inverse enthalpy function given the internal energy, latent heat content, and heat
capacity under the free water freezing characteristic.
"""
@inline function energy_to_temperature(::FreeWater, U::NF, Lθ::NF, C::NF) where {NF}
    return if U < -Lθ
        # Case 1: U < -Lθ → frozen
        (U + Lθ) / C
    elseif U >= zero(U)
        # Case 2a: U ≥ 0 → thawed
        U / C
    else
        # Case 2b: -Lθ ≤ U < 0 → phase change
        zero(NF)
    end
    # One-liner version:
    # return (U < -Lθ)*(U + Lθ) / C
end
