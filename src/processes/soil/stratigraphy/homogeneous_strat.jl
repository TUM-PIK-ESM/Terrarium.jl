"""
    $TYPEDEF

Represents a soil stratigraphy of well mixed material with homogeneous soil texture.

Properties:
$TYPEDFIELDS
"""
struct HomogeneousStratigraphy{NF, SoilPorosity<:AbstractSoilPorosity{NF}} <: AbstractStratigraphy{NF}
    "Material composition of mineral soil component"
    texture::SoilTexture{NF}
    
    "Parameterization of soil porosity"
    porosity::SoilPorosity
end

function HomogeneousStratigraphy(
    ::Type{NF};
    texture::AbstractSoilTexture{NF} = SoilTexture(NF),
    porosity::AbstractSoilPorosity{NF} = ConstantSoilPorosity(NF)
) where {NF}
    return HomogeneousStratigraphy(texture, porosity)
end

# Kernel functions

soil_texture(i, j, k, grid, state, strat::HomogeneousStratigraphy) = strat.texture

@inline function organic_fraction(
    i, j, k, grid, state,
    strat::HomogeneousStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    ρ_soc = density_soc(i, j, k, grid, state, bgc)
    ρ_org = density_pure_soc(bgc)
    por = organic_porosity(i, j, k, grid, state, strat.porosity, strat.texture)
    organic = ρ_soc / ((1 - por) * ρ_org)
    return organic
end

"""
    $TYPEDSIGNATURES

Compute the porosity of the soil volume at the given indices.
"""
@propagate_inbounds function porosity(
    i, j, k, grid, state,
    strat::HomogeneousStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    organic = organic_fraction(i, j, k, grid, state, bgc, strat)
    texture = strat.soil_texture
    por_m = mineral_porosity(i, j, k, grid, state, strat.porosity, texture)
    por_o = organic_porosity(i, j, k, grid, state, strat.porosity, texture)
    return (1 - organic)*por_m + organic*por_o
end

"""
    $SIGNATURES

Construct a `SoilVolume` object summarizing the material composition of the soil volume
at the given indices `i, j, k` on `grid`.
"""
@propagate_inbounds function soil_volume(
    i, j, k, grid, state,
    strat::AbstractStratigraphy,
    hydrology::AbstractSoilHydrology,
    bgc::AbstractSoilBiogeochemistry
)
    # get current saturation state and liquid fraction
    sat = saturation_water_ice(i, j, k, grid, state, hydrology)
    liq = liquid_water_fraction(i, j, k, grid, state, hydrology)
    por = porosity(i, j, k, grid, state, strat, bgc)
    solid = soil_matrix(i, j, k, grid, state, strat, bgc)
    return SoilVolume(por, sat, liq, solid)
end

"""
    $TYPEDSIGNATURES

Compute and return the soil solid matrix at index `i, j, k` on `grid`. The default
implementation assumes a simple `MineralOrganic` parameterization of the solid material.
"""
@propagate_inbounds function soil_matrix(
    i, j, k, grid, state,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    organic = organic_fraction(i, j, k, grid, state, strat, bgc)
    texture = soil_texture(i, j, k, grid, state, strat)
    return MineralOrganic(texture, organic)
end
