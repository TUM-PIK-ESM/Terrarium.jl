@kwdef struct AutotrophicRespirationModel{NF} <: AbstractAutotrophicRespiration
    # TODO add autotrophic respiration parameters

end

variables(autoresp::AutotrophicRespirationModel) = (
    auxiliary(:Ra, XY()), # Autotrophic respiration (Ra) kgC/m²/day
    #TODO should NPP be computed here?
    auxiliary(:NPP, XY()), # Net Primary Production (NPP) kgC/m²/day
)

function compute_auxiliary!(idx, state, model, autoresp::AutotrophicRespirationModel)
    i, j = idx

    # Compute Ra
    # Placeholder assuming a simple linear relationship with C_veg
    state.Ra[i, j] = 0.05 * state.veg_carbon_pool[i, j] 

    # Compute NPP
    state.NPP[i, j] = state.GPP[i, j] - state.Ra[i, j]
end

   