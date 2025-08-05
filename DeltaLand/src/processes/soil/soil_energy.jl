abstract type AbstractHeatOperator end
struct ExplicitHeatConduction <: AbstractHeatOperator end

@kwdef struct SoilEnergyBalance{
    HeatOperator<:AbstractHeatOperator,
    ThermalProps<:SoilThermalProperties,
} <: AbstractSoilEnergyBalance
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

@inline function thermalconductivity(idx, state, model::AbstractSoilModel, energy::SoilEnergyBalance)
    conds = getproperties(energy.thermal_properties.cond)
    fracs = soil_volumetric_fractions(idx, state, model)
    # apply bulk conductivity weighting
    return energy.thermal_properties.cond_bulk(conds, fracs)
end

@inline function heatcapacity(idx, state, model::AbstractSoilModel, energy::SoilEnergyBalance)
    heatcaps = getproperties(energy.thermal_properties.heatcap)
    fracs = soil_volumetric_fractions(idx, state, model)
    # for heat capacity, we just do a weighted average
    return sum(fastmap(*, heatcaps, fracs))
end

@inline function compute_auxiliary!(idx, state, model, energy::SoilEnergyBalance)
    # currently nothing to do here
    return nothing
end

@inline function compute_tendencies!(idx, state, model, energy::SoilEnergyBalance)
    i, j, k = idx

    # Get underlying Oceananigans grid
    grid = get_field_grid(get_grid(model))

    # Get temperature state
    T = state.temperature

    # Get thermal conductivity
    κ = thermalconductivity(idx, state, model, energy)
    
    # Interior heat fluxes
    dUdt = -∂zᵃᵃᶜ(i, j, k, grid, diffusive_flux, κ, T)
    
    # Ground heat flux (upper boundary)
    Qg = @inbounds state.ground_heat_flux[i, j, k]
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    dUdt += ifelse(k == grid.Nz, -Qg / Δz, zero(grid))

    # Geothermal heat flux (lower boundary)
    # TODO: should this just be a flux boundary condition on the tendency?
    Qgeo = @inbounds state.geothermal_heat_flux[i, j, k]
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    dUdt += ifelse(k == 1, Qgeo / Δz, zero(grid))

    # Accumulate tendency
    state.internal_energy_tendency[i, j, k] += dUdt
end

# Diffusive heat flux term passed to ∂z operator
@inline diffusive_flux(i, j, k, grid, κ, T) = -κ * ∂zᵃᵃᶠ(i, j, k, grid, T)

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

function closure!(state, model::AbstractSoilModel, closure::TemperatureEnergyClosure)
    grid = get_field_grid(get_grid(model))
    launch!(grid, _closure_kernel, state, model, closure)
    return nothing
end

function invclosure!(state, model::AbstractSoilModel, ::TemperatureEnergyClosure)
    grid = get_grid(model)
    launch!(grid, _invclosure_kernel, state, model, closure)
    return nothing
end

@kernel function _closure_kernel(state, model::AbstractSoilModel, ::TemperatureEnergyClosure)
    idx = @index(Global, NTuple)
    temperature_to_energy!(idx, state, model, get_freezecurve(get_soil_hydrology(model)))
end

@kernel function _invclosure_kernel(state, model, ::TemperatureEnergyClosure)
    idx = @index(Global, NTuple)
    energy_to_temperature!(idx, state, model, get_freezecurve(get_soil_hydrology(model)))
end

@inline function temperature_to_energy!(idx, state, model, ::FreezeCurves.FreeWater)
    i, j, k = idx
    constants = get_constants(model)
    T = state.temperature[i, j, k] # assumed given
    C = heatcapacity(idx, state, model, get_soil_energy_balance(model))
    L = constants.ρw*constants.Lsl
    por = porosity(idx, state, model, get_soil_hydrology(model))
    sat = state.pore_water_ice_saturation[i, j, k]
    θwi = sat*por
    Lθ = L*θwi
    # calculate unfrozen water content from temperature
    # N.B. For the free water freeze curve, the mapping from temperature to unfrozen water content
    # within the phase change region is indeterminate since it is assumed that T = 0. As such, we
    # have to assume here that the liquid water fraction is zero if T <= 0. This method shold therefore
    # only be used for initialization and should **not** be involved in the calculation of tendencies.
    liquid_water_frac = ifelse(
        T > zero(T),
        θwi,
        zero(θwi),
    )
    state.liquid_water_fraction[i, j, k] = liquid_water_frac
    # compute energy from temperature, heat capacity, and ice fraction
    U = state.internal_energy[i, j, k] = T*C - L*θwi*(1 - liquid_water_frac)
    return U
end

@inline function energy_to_temperature!(idx, state, model, ::FreezeCurves.FreeWater)
    i, j, k = idx
    constants = get_constants(model)
    U = state.internal_energy[i, j, k] # assumed given
    C = heatcapacity(idx, state, model, get_soil_energy_balance(model))
    L = constants.ρw*constants.Lsl
    por = porosity(idx, state, model, get_soil_hydrology(model))
    sat = state.pore_water_ice_saturation[i, j, k]
    θwi = sat*por
    Lθ = L*θwi
    # calculate unfrozen water content
    liquid_water_frac, _ = FreezeCurves.freewater(U, one(θwi), L)
    state.liquid_water_fraction[i, j, k] = liquid_water_frac
    # calculate temperature from internal energy and liquid water fraction
    T = state.temperature[i, j, k] = ifelse(
        U < -Lθ,
        # Case 1: U < 0 → frozen
        (U - Lθ) / C,
        # Case 2: U ≥ -Lθ
        ifelse(
            U >= zero(U),
            # Case 2a: U ≥ 0 → thawed
            U / C,
            # Case 2b: -Lθ ≤ U < 0 → phase change
            zero(eltype(state.temperature))
        )
    )
    return T
end
