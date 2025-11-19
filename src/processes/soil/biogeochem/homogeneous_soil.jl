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

soil_texture(i, j, k, state, strat::HomogeneousSoil) = strat.texture

variables(strat::HomogeneousSoil) = ()

# do nothing for now
@inline initialize!(state, model, strat::HomogeneousSoil) = nothing

@inline compute_auxiliary!(state, model, strat::HomogeneousSoil) = nothing

@inline compute_tendencies!(state, model, strat::HomogeneousSoil) = nothing

"""
Compute and return a `SoilComposition` summarizing the material composition of the soil volume at the given indices.
"""
@inline function soil_composition(
        i, j, k, state,
        strat::HomogeneousSoil,
        hydrology::AbstractSoilHydrology,
        bgc::AbstractSoilBiogeochemistry
    )
    # get current saturation state and liquid fraction
    sat = state.saturation_water_ice[i, j, k]
    liq = state.liquid_water_fraction[i, j, k]
    por = porosity(i, j, k, state, hydrology, strat, bgc)
    # there is some slight redundant computation here; consider merging into one method?
    org = organic_fraction(i, j, k, state, bgc)
    return SoilComposition(por, sat, liq, org, strat.texture)
end
