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

variables(strat::HomogeneousSoil) = (
    auxiliary(:porosity, XYZ(), domain=UnitInterval(), desc="Bulk porosity (void fraction) of the soil volume"),
)

@inline initialize!(state, model, strat::HomogeneousSoil) = nothing

function compute_auxiliary!(state, model, strat::HomogeneousSoil)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    bgc = get_soil_biogeochemistry(model)
    launch!(state, grid, :xyz, compute_porosity!, strat, hydrology, bgc)
end

@inline compute_tendencies!(state, model, strat::HomogeneousSoil) = nothing

# Kernels

@kernel function compute_porosity!(state, grid, strat::HomogeneousSoil, hydrology::AbstractSoilHydrology, bgc::AbstractSoilBiogeochemistry)
    i, j, k = @index(Global, NTuple)
    state.porosity[i, j, k] = compute_porosity(i, j, k, state, strat, hydrology, bgc)
end

# Kernel functions

porosity(i, j, k, state, ::AbstractStratigraphy) = @inbounds state.porosity[i, j, k]

"""
    compute_porosity(i, j, k, state, strat::HomogeneousSoil, bgc::AbstractSoilBiogeochemistry)

Compute the porosity of the soil volume at the given indices based on the current state as well as the stratigraphy, and biogeochemistry configurations.
"""
@inline function compute_porosity(i, j, k, state, strat::HomogeneousSoil, hydrology::AbstractSoilHydrology, bgc::AbstractSoilBiogeochemistry)
    org = organic_fraction(i, j, k, state, bgc)
    texture = soil_texture(i, j, k, state, strat)
    props = get_hydraulic_properties(hydrology)
    return (1 - org)*mineral_porosity(props, texture) + org*organic_porosity(i, j, k, state, bgc)
end

"""
Compute and return a `SoilComposition` object summarizing the material composition of the soil volume at the given indices.
"""
@inline function soil_composition(
    i, j, k, state,
    strat::HomogeneousSoil,
    energy::AbstractSoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    bgc::AbstractSoilBiogeochemistry
)
    # get current saturation state and liquid fraction
    sat = saturation_water_ice(i, j, k, state, hydrology)
    liq = liquid_water_fraction(i, j, k, state, energy)
    por = porosity(i, j, k, state, strat)
    # there is some slight redundant computation here; consider merging into one method?
    org = organic_fraction(i, j, k, state, bgc)
    return SoilComposition(por, sat, liq, org, strat.texture)
end
