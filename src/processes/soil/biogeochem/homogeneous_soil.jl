"""
    $TYPEDEF

Represents a soil stratigraphy of well mixed material with homogeneous soil texture.

Properties:
$TYPEDFIELDS
"""
struct HomogeneousSoil{NF} <: AbstractStratigraphy{NF}
    "Material composition of mineral soil componnet"
    texture::SoilTexture{NF}
end

HomogeneousSoil(::Type{NF}; texture::SoilTexture{NF} = SoilTexture(NF, :sand)) where {NF} = HomogeneousSoil{NF}(texture)

soil_texture(idx, state, strat::HomogeneousSoil) = strat.texture

variables(strat::HomogeneousSoil) = ()

# do nothing for now
@inline initialize!(state, model, strat::HomogeneousSoil) = nothing

@inline compute_auxiliary!(state, model, strat::HomogeneousSoil) = nothing

@inline compute_tendencies!(state, model, strat::HomogeneousSoil) = nothing

"""
Compute and return a `SoilComposition` summarizing the material composition of the soil volume at index `idx`.
"""
@inline function soil_composition(
    idx, state,
    strat::HomogeneousSoil,
    hydrology::AbstractSoilHydrology,
    bgc::AbstractSoilBiogeochemistry
)
    # get current saturation state and liquid fraction
    sat = state.saturation_water_ice[idx...]
    liq = state.liquid_water_fraction[idx...]
    por = porosity(idx, state, hydrology, strat, bgc)
    # there is some slight redundant computation here; consider merging into one method?
    org = organic_fraction(idx, state, bgc)
    return SoilComposition(por, sat, liq, org, strat.texture)
end
