"""
    $TYPEDEF
"""
@kwdef struct HomogeneousSoil{NF} <: AbstractStratigraphy
    "Material composition of mineral soil componnet"
    texture::SoilTexture{NF,NF,NF} = SoilTexture(:sand)
end

get_soil_texture(strat::HomogeneousSoil) = strat.texture

variables(strat::HomogeneousSoil) = ()

# do nothing for now
@inline update_state!(idx, state, model, strat::HomogeneousSoil) = nothing

@inline compute_tendencies!(idx, state, model, strat::HomogeneousSoil) = nothing
