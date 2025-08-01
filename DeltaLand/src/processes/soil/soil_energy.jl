abstract type AbstractHeatOperator end
struct VerticalHeatConduction <: AbstractHeatOperator end

@kwdef struct SoilEnergyBalance{
    HeatOperator<:AbstractHeatOperator,
    ThermalProps<:SoilThermalProperties,
} <: AbstractSoilEnergyBalance
    "Heat transport operator"
    operator::HeatOperator = VerticalHeatConduction()

    # Note: Would it make more sense for these properties to be defined in the stratigraphy?
    "Soil thermal properties"
    thermal_properties::ThermalProps = SoilThermalProperties()
end

variables(::SoilEnergyBalance) = (
    prognostic(:temperature, XYZ(), closure=TemperatureEnergyClosure()),
    auxiliary(:pore_water_ice_saturation, XYZ()),
    axuiliary(:liquid_water_fraction, XYZ()),
    auxiliary(:ground_heat_flux, XY()),
    auxiliary(:geothermal_heat_flux, XY()),
)

@inline function thermalconductivity(i, j, k, state, model::AbstractSoilModel, energy::SoilEnergyBalance)
    conds = getproperties(energy.thermal_properties.cond)
    fracs = soil_volumetric_fractions(i, j, k, state, model)
    # apply bulk conductivity weighting
    return energy.thermal_properties.cond_bulk(conds, fracs)
end

@inline function heatcapacity(i, j, k, state, model::AbstractSoilModel, energy::SoilEnergyBalance)
    heatcaps = getproperties(energy.thermal_properties.heatcap)
    fracs = soil_volumetric_fractions(i, j, k, state, model)
    # for heat capacity, we just do a weighted average
    return sum(fastmap(*, heatcaps, fracs))
end

@inline function update_state(i, j, k, state, model, energy::SoilEnergyBalance)
    enthalpyinv(i, j, k, state, model, freezecurve(model))
end

@inline function compute_tendencies(i, j, k, state, model, energy::SoilEnergyBalance)
    # Get underlying Oceananigans grid
    grid = get_field_grid(get_grid(model))

    # Get temperature state
    T = state.temperature

    # Get thermal conductivity
    κ = thermalconductivity(i, j, k, state, model, energy)
    
    # Interior heat flux
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

@inline diffusive_flux(i, j, k, grid, κ, T) = -κ * ∂zᵃᵃᶠ(i, j, k, grid, T)

"""
    TemperatureEnergyClosure

Defines the constitutive relationship between temperature and the internal energy of the system, i.e:

```math
h(T) = T\times C(T) + L_f \theta(T)
```
where 
"""
struct TemperatureEnergyClosure end

varname(::TemperatureEnergyClosure) = :internal_energy

@inline function enthalpyinv(i, j, k, state, model, ::FreezeCurves.FreeWater)
    constants = get_constants(model)
    let H = state.enthalpy[i, j, k], # assumed prognostic
        C = heatcapacity(i, j, k, state, model, get_energy_balance(model)),
        L = constants.ρw*constants.Lsl,
        por = porosity(i, j, k, state, model, get_stratigraphy(model)),
        sat = state.pore_water_ice_saturation[i, j, k],
        θwi = sat*por,
        Lθ = L*θwi;
        # calculate unfrozen water content
        liquid_water_frac, _ = FreezeCurves.freewater(H, one(θwi), L)
        state.liquid_water_fraction[i, j, k] = liquid_water_frac
        # calculate temperature
        T = ifelse(
            H < zero(θwi),
            # Case 1: H < 0 -> frozen
            H / C,
            # Case 2: H >= 0
            ifelse(
                H >= Lθ,
                # Case 2a: H >= Lθ -> thawed
                (H - Lθ) / C,
                # Case 2b: 0 <= H < Lθ -> phase change
                zero(eltype(state.temperature))
            )
        )
        return T
    end
end
