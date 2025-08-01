struct ImmobileSoilWater <: AbstractSoilHydrology end

variables(::ImmobileSoilWater) = ()

@inline function update_state(i, j, k, state, model::AbstractSoilModel, hydrology::ImmobileSoilWater)
    strat = get_stratigraphy(model)
    # set saturation level of pore water/ice to value specified by stratigraphy
    state.pore_water_ice_saturation[i,j,k] = pore_water_ice_saturation(i, j, k, state, model, strat)
end

@inline compute_tendencies(i, j, k, state, model, strat::ImmobileSoilWater) = nothing
