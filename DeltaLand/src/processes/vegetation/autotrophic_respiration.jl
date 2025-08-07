# TODO maybe change the name later, if the PALADYN autotrophic respiration approach has a more specific name

@kwdef struct PALADYNAutotrophicRespiration{NF} <: AbstractAutotrophicRespiration
    # TODO add autotrophic respiration parameters
end

variables(autoresp::PALADYNAutotrophicRespiration) = (
    auxiliary(:GPP, XY()), # Gross Primary Production [kgC/m²/day]
    auxiliary(:C_veg, XY()), # Vegetation carbon pool [kgC/m²]
    auxiliary(:Ra, XY()), # Autotrophic respiration [kgC/m²/day]
    auxiliary(:NPP, XY()), # Net Primary Production [kgC/m²/day]
)

@inline function compute_auxiliary!(idx, state, model::AbstractVegetationModel, autoresp::PALADYNAutotrophicRespiration{NF}) where NF
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

   