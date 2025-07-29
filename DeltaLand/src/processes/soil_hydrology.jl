struct ImmobileSoilWater <: AbstractSoilHydrology end

@inline update_state(i, j, k, state, model, hydrology::ImmobileSoilWater) = nothing
