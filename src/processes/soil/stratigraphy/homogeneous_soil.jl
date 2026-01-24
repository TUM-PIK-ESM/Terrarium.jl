"""
    $TYPEDEF

Represents a soil stratigraphy of well mixed material with homogeneous soil texture.

Properties:
$TYPEDFIELDS
"""
struct HomogeneousSoil{NF} <: AbstractStratigraphy{NF}
    "Material composition of mineral soil component"
    texture::SoilTexture{NF}
end

HomogeneousSoil(::Type{NF}; texture::SoilTexture{NF} = SoilTexture(NF, :sand)) where {NF} = HomogeneousSoil{NF}(texture)

soil_texture(i, j, k, grid, state, strat::HomogeneousSoil) = strat.texture

variables(strat::HomogeneousSoil) = (
    auxiliary(:porosity, XYZ(), domain=UnitInterval(), desc="Bulk porosity (void fraction) of the soil volume"),
)

@inline function initialize!(state, model, strat::HomogeneousSoil)
    compute_auxiliary!(state, model, strat)
end

function compute_auxiliary!(state, model, strat::HomogeneousSoil)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    bgc = get_soil_biogeochemistry(model)
    launch!(grid, XYZ, compute_porosity!, state, strat, hydrology, bgc)
end

@inline compute_tendencies!(state, model, strat::HomogeneousSoil) = nothing

# Kernels

@kernel function compute_porosity!(state, grid, strat::HomogeneousSoil, hydrology::AbstractSoilHydrology, bgc::AbstractSoilBiogeochemistry)
    i, j, k = @index(Global, NTuple)
    state.porosity[i, j, k] = compute_porosity(i, j, k, grid, state, strat, hydrology, bgc)
end

# Kernel functions

porosity(i, j, k, grid, state, ::AbstractStratigraphy) = @inbounds state.porosity[i, j, k]

"""
    compute_porosity(i, j, k, grid, state, strat::HomogeneousSoil, bgc::AbstractSoilBiogeochemistry)

Compute the porosity of the soil volume at the given indices based on the current state as well as the stratigraphy, and biogeochemistry configurations.
"""
@inline function compute_porosity(i, j, k, grid, state, strat::HomogeneousSoil, hydrology::AbstractSoilHydrology, bgc::AbstractSoilBiogeochemistry)
    org = organic_fraction(i, j, k, grid, state, bgc)
    texture = soil_texture(i, j, k, grid, state, strat)
    props = get_hydraulic_properties(hydrology)
    return (1 - org)*mineral_porosity(props, texture) + org*organic_porosity(i, j, k, grid, state, bgc)
end

"""
    $SIGNATURES

Construct a `SoilVolume` object summarizing the material composition of the soil volume.
"""
@inline function soil_volume(
    i, j, k, grid, state,
    strat::HomogeneousSoil,
    hydrology::AbstractSoilHydrology,
    bgc::AbstractSoilBiogeochemistry
)
    # get current saturation state and liquid fraction
    sat = saturation_water_ice(i, j, k, grid, state, hydrology)
    liq = liquid_water_fraction(i, j, k, grid, state, hydrology)
    por = porosity(i, j, k, grid, state, strat)
    solid = @inbounds soil_matrix(i, j, k, grid, state, strat, bgc)
    return SoilVolume(por, sat, liq, solid)
end

# Soil matrix for homogeneous soil stratigraphy
@inline function soil_matrix(
    i, j, k, grid, state,
    strat::HomogeneousSoil,
    bgc::AbstractSoilBiogeochemistry
)
    org = organic_fraction(i, j, k, grid, state, bgc)
    tex = strat.texture
    return MineralOrganic(tex, org)
end
