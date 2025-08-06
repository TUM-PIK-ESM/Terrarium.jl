@kwdef struct LUEPhotosynthesis{NF} <: AbstractPhotosynthesis
    # TODO add photosynthesis parameters 
end

variables(photo::LUEPhotosynthesis) = (
    auxiliary(:GPP, XY()), # Gross Primary Production (GPP) kgC/m²/day
)

@inline function compute_auxiliary!(idx, state, model::AbstractVegetationModel, photo::LUEPhotosynthesis)
    i, j = idx

    # TODO add photosynthesis implementation
    # Needs state.λc from stomatal conductance
    
    # For now, set GPP to a random value 
    state.GPP[i, j] = 0.01*rand() 
end
