abstract type AbstractHeatOperator end
struct ExplicitHeatConduction <: AbstractHeatOperator end

"""
    $TYPEDEF

Standard implementation of the soil energy balance accounting for freezing and thawing of pore water/ice.

Properties:
$TYPEDFIELDS
"""
@kwdef struct SoilEnergyBalance{
    NF,
    HeatOperator<:AbstractHeatOperator,
    ThermalProps<:SoilThermalProperties{NF},
} <: AbstractSoilEnergyBalance{NF}
    "Heat transport operator"
    operator::HeatOperator = ExplicitHeatConduction()

    "Soil thermal properties"
    thermal_properties::ThermalProps = SoilThermalProperties()
end

variables(::SoilEnergyBalance) = (
    prognostic(:temperature, XYZ(), TemperatureEnergyClosure()),
    auxiliary(:pore_water_ice_saturation, XYZ()),
    auxiliary(:liquid_water_fraction, XYZ()),
    auxiliary(:ground_heat_flux, XY()),
    auxiliary(:geothermal_heat_flux, XY()),
)

compute_auxiliary!(state, model, energy::SoilEnergyBalance) = nothing

function compute_tendencies!(state, model, energy::SoilEnergyBalance)
    grid = get_grid(model)
    field_grid = get_field_grid(grid)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    launch!(grid, :xyz, compute_energy_tendency!, state, field_grid, energy, hydrology, strat, bgc)
    return nothing
end

@kernel function compute_energy_tendency!(
    state,
    grid,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    idx = @index(Global, NTuple)
    state.internal_energy_tendency[idx...] += energy_tendency(idx, state, grid, energy, hydrology, strat, bgc)
end

@inline function energy_tendency(
    idx, state, grid,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    i, j, k = idx

    # Get temperature field
    T = state.temperature
    
    # Retrieve ground heat flux (upper boundary)
    # Negative since fluxes are always positive upwards
    Qg = @inbounds state.ground_heat_flux[i, j, k]
    # Retrieve geothermal heat flux (lower boundary)
    Qgeo = @inbounds state.geothermal_heat_flux[i, j, k]

    return (
        # Interior heat fluxes
        -∂zᵃᵃᶜ(i, j, k, grid, diffusive_heat_flux, state, energy, hydrology, strat, bgc)
        # Upper boundary flux
        + energy_boundary_tendency(i, j, k, grid, grid.Nz, -Qg)
        # Lower boundary flux
        + energy_boundary_tendency(i, j, k, grid, 1, Qgeo)
    )
end

@inline function thermalconductivity(
    i, j, k, grid, state,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
)
    fracs = soil_volumetric_fractions((i, j, k), state, strat, hydrology, bgc)
    return thermalconductivity(energy.thermal_properties, fracs)
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

# Tendency for boundary cells; zero for all interior grid cells
@inline function energy_boundary_tendency(i, j, k, grid, I, Q)
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    return (k == I) * (Q / Δz)
end

"""
    TemperatureEnergyClosure

Defines the constitutive relationship between temperature and the internal energy, U, of the system, i.e:

```math
U(T) = T\\times C(T) - L_f \\theta_{wi} (1 - F(T))
```
where T is temperature, C(T) is the temperature-dependent heat capacity, L_f is the
volumetric latent heat of fusion, and F(T) is the constitutive relation between temperature
and the unfrozen fraction of pore water. Note that this formulation implies that the zero
"""
struct TemperatureEnergyClosure <: AbstractClosureRelation end

varname(::TemperatureEnergyClosure) = :internal_energy

function closure!(state, model::AbstractSoilModel, ::TemperatureEnergyClosure)
    grid = get_grid(model)
    energy = get_soil_energy_balance(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    constants = get_constants(model)
    launch!(grid, :xyz, temperature_to_energy!, state, energy, hydrology, strat, bgc, constants)
    return nothing
end

function invclosure!(state, model::AbstractSoilModel, ::TemperatureEnergyClosure)
    grid = get_grid(model)
    energy = get_soil_energy_balance(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    constants = get_constants(model)
    launch!(grid, :xyz, energy_to_temperature!, state, energy, hydrology, strat, bgc, constants)
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
    idx = @index(Global, NTuple)
    fc = get_freezecurve(hydrology)
    temperature_to_energy!(idx, state, fc, energy, hydrology, strat, bgc, constants)
end

@inline function temperature_to_energy!(
    idx, state, ::FreezeCurves.FreeWater,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
)
    i, j, k = idx
    T = state.temperature[i, j, k] # assumed given
    L = constants.ρw*constants.Lsl
    por = porosity(idx, state, hydrology, strat, bgc)
    sat = state.pore_water_ice_saturation[i, j, k]
    Lθ = L*sat*por
    # calculate unfrozen water content from temperature
    # N.B. For the free water freeze curve, the mapping from temperature to unfrozen water content
    # within the phase change region is indeterminate since it is assumed that T = 0. As such, we
    # have to assume here that the liquid water fraction is zero if T <= 0. This method should therefore
    # only be used for initialization and should **not** be involved in the calculation of tendencies.
    liq = state.liquid_water_fraction[i, j, k] = ifelse(
        T > zero(T),
        sat,
        zero(sat),
    )
    fracs = soil_volumetric_fractions(idx, state, strat, hydrology, bgc)
    C = heatcapacity(energy.thermal_properties, fracs)
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
    idx = @index(Global, NTuple)
    fc = get_freezecurve(hydrology)
    energy_to_temperature!(idx, state, fc, energy, hydrology, strat, bgc, constants)
end

@inline function energy_to_temperature!(
    idx, state, ::FreezeCurves.FreeWater,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
)
    i, j, k = idx
    U = state.internal_energy[i, j, k] # assumed given
    L = constants.ρw*constants.Lsl
    por = porosity(idx, state, hydrology, strat, bgc)
    sat = state.pore_water_ice_saturation[i, j, k]
    Lθ = L*sat*por
    # calculate unfrozen water content
    state.liquid_water_fraction[i, j, k] = clamp(one(sat) + U/Lθ, zero(sat), one(sat))

    fracs = soil_volumetric_fractions(idx, state, strat, hydrology, bgc)
    C = heatcapacity(energy.thermal_properties, fracs)
    # calculate temperature from internal energy and liquid water fraction
    # Case 1: U < -L0 -> frozen; Case 2: U ≥ -L0; Case 2a: U ≥ 0 -> thawed; Case 2b: -L0 ≤ U < 0 -> phase change                            
    T = state.temperature[i, j, k] = (U < -Lθ) * (U + Lθ) / C + (U >= zero(U)) * U / C
    return T
end
