"""
    $TYPEDEF

Represents a soil stratigraphy of well mixed material with homogeneous soil texture.

Properties:
$TYPEDFIELDS
"""
@kwdef struct HomogeneousSoil{NF} <: AbstractStratigraphy{NF}
    "Material composition of mineral soil componnet"
    texture::SoilTexture{NF,NF,NF} = SoilTexture(:sand)
end

get_soil_texture(strat::HomogeneousSoil) = strat.texture

variables(strat::HomogeneousSoil) = ()

# do nothing for now
@inline compute_auxiliary!(state, model, strat::HomogeneousSoil) = nothing

@inline compute_tendencies!(state, model, strat::HomogeneousSoil) = nothing
