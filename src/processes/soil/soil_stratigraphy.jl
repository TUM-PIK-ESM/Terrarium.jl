"""
    $TYPEDEF

Represents a soil stratigraphy of well mixed material with homogeneous soil texture.

Properties:
$TYPEDFIELDS
"""
@kwdef struct HomogeneousSoil{NF} <: AbstractStratigraphy{NF}
    "Material composition of mineral soil componnet"
    texture::SoilTexture{NF} = SoilTexture(:sand)
end

get_soil_texture(strat::HomogeneousSoil) = strat.texture

variables(strat::HomogeneousSoil) = ()

# do nothing for now
@inline compute_auxiliary!(state, model, strat::HomogeneousSoil) = nothing

@inline compute_tendencies!(state, model, strat::HomogeneousSoil) = nothing

@inline function soil_characteristic_fractions(
    idx, state,
    strat::AbstractStratigraphy,
    hydrology::AbstractSoilHydrology,
    bgc::AbstractSoilBiogeochemistry
)
    sat = state.pore_water_ice_saturation[idx...]
    por = porosity(idx, state, hydrology, strat, bgc)
    ## there is some slight redundant computation here; consider merging into one method?
    org = organic_fraction(idx, state, bgc)
    return (; sat, por, org)
end

@inline function soil_volumetric_fractions(
    idx, state,
    strat::AbstractStratigraphy,
    hydrology::AbstractSoilHydrology,
    bgc::AbstractSoilBiogeochemistry
)
    # get characteristic fractions
    (; sat, por, org) = soil_characteristic_fractions(idx, state, strat, hydrology, bgc)
    # get fraction of unfrozen pore water
    liq = state.liquid_water_fraction[idx...]
    # calculate volumetric fractions
    water_ice = sat*por
    water = water_ice*liq
    ice = water_ice*(1-liq)
    air = (1-sat)*por
    mineral = (1-por)*(1-org)
    organic = (1-por)*org
    return (; water, ice, air, mineral, organic)
end
