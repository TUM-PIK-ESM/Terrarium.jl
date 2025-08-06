# TODO maybe change the name later, if the PALADYN 
# autotrophic respiration approach has a more specific name

@kwdef struct PaladynAutotrophicRespiration{NF} <: AbstractAutotrophicRespiration
    # TODO add autotrophic respiration parameters
end

variables(autoresp::PaladynAutotrophicRespiration) = (
    auxiliary(:GPP, XY()), # Gross Primary Production (GPP) kgC/m²/day
    auxiliary(:C_veg, XY()), # Vegetation carbon pool (C_veg) kgC/m²
    auxiliary(:Ra, XY()), # Autotrophic respiration (Ra) kgC/m²/day
    auxiliary(:NPP, XY()), # Net Primary Production (NPP) kgC/m²/day
)

@inline function compute_auxiliary!(idx, state, model::AbstractVegetationModel, autoresp::PaladynAutotrophicRespiration{NF}) where NF
    i, j = idx

    # Compute maintenance respiration Rm
    # TODO add maintenance respiration implementaion
    # Placeholder assuming a simple linear relationship with C_veg
    Rm = NF(0.05) * state.C_veg[i, j] 

    # Compute growth respiration Rg
    Rg = NF(0.25) * (state.GPP[i, j] - Rm)

    # Compute autotrophic respiration Ra
    state.Ra[i, j] = Rm + Rg

    # Compute NPP
    state.NPP[i, j] = state.GPP[i, j] - state.Ra[i, j]
end

   