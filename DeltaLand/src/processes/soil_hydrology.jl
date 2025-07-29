struct ImmobileSoilWater <: AbstractSoilHydrology end

@inline update_state(i, j, k, state, model, hydrology::ImmobileSoilWater) = nothing

@inline compute_tendencies(i, j, k, state, model, strat::ImmobileSoilWater) = nothing
