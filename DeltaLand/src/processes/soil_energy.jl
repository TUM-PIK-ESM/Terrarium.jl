abstract type AbstractHeatOperator end
struct VerticalHeatConduction <: AbstractHeatOperator end

# TODO: In principle, these types could change for different soil parameterizations.
# This is something we should ideally allow for.
@kwdef struct SoilThermalProperties{NF,CondWeighting}
    "Thermal conductivities for all constituents"
    cond::SoilThermalConductivities{NF} = SoilThermalConductivities()

    "Thermal conductivity mixing scheme"
    cond_bulk::CondWeighting = InverseQuadratic()

    "Thermal conductivities for all constituents"
    heatcap::SoilHeatCapacities{NF} = SoilHeatCapacities()
end

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

function initialize(::SoilEnergyBalance, grid::AbstractLandGrid, initializer)
    # TODO: it would be cool if we could formally specify constitutive relations like temperature <--> enthalpy
    enthalpy = initialize3D(initializer, grid, :enthalpy)
    enthalpy_tendency = initialize3D(initializer, grid, :enthalpy_tendency)
    temperature = initialize3D(initializer, grid, :temperature)
    pore_water_ice_saturation = initialize3D(initializer, grid, :pore_water_ice_saturation)
    liquid_water_fraction = initialize3D(initializer, grid, :liquid_water_fraction)
    # 2D boundary variables
    ground_heat_flux = initialize2D(initializer, grid, :ground_heat_flux)
    geothermal_heat_flux = initialize2D(initializer, grid, :geothermal_heat_flux)
    # return state named tuples
    prognostic = (; enthalpy)
    tendencies = (; enthalpy_tendency)
    auxiliary = (; temperature, pore_water_ice_saturation, liquid_water_fraction)
    boundary_vars = (; ground_heat_flux, geothermal_heat_flux)
    return StateVariables(; prognostic, tendencies, auxiliary, boundary_vars)
end

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
        state.temperature[i, j, k] = ifelse(
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
    end
end

@inline function update_state(i, j, k, state, model, energy::SoilEnergyBalance)
    enthalpyinv(i, j, k, state, model, freezecurve(model))
end

@inline function compute_tendencies(i, j, k, state, model, energy::SoilEnergyBalance)
    # Get underlying Oceananigans grid
    grid = get_grid_impl(get_grid(model))

    # Get temperature state
    T = state.temperature

    # Get thermal conductivity
    κ = thermalconductivity(i, j, k, state, model, energy)
    
    # Interior heat flux
    dHdt = -∂zᵃᵃᶜ(i, j, k, grid, diffusive_flux, κ, T)
    
    # Ground heat flux (upper boundary)
    Qg = @inbounds state.ground_heat_flux[i, j, k]
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    dHdt += ifelse(k == grid.Nz, -Qg / Δz, zero(grid))

    # Geothermal heat flux (lower boundary)
    # TODO: should this just be a flux boundary condition on the tendency?
    Qgeo = @inbounds state.geothermal_heat_flux[i, j, k]
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    dHdt += ifelse(k == 1, Qgeo / Δz, zero(grid))

    # Accumulate tendencies
    state.enthalpy_tendency .+= dHdt
end

@inline diffusive_flux(i, j, k, grid, κ, T) = -κ * ∂zᵃᵃᶠ(i, j, k, grid, T)
