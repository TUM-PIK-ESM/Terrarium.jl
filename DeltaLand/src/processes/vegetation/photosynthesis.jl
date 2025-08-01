@kwdef struct PhotosynthesisFarquhar{NF} <: AbstractPhotosynthesis
    # parameters here
end

variables(photo::PhotosynthesisFarquhar) = (
    auxiliary(:GPP, XY()),
)

function update_state!(idx, state, model, photo::PhotosynthesisFarquhar)
    i, j, k = idx

    # Add photosynthesis impl here

    # state.GPP[i, j, k] = ...
end
