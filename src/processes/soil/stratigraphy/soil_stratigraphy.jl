"""
    $TYPEDEF

Represents a soil stratigraphy as a stack of named [soil horizons](https://en.wikipedia.org/wiki/Soil_horizon).
Each soil horizon is assumed to have internally homogeneous soil properties. The number of horizons and their
respective names are defined by the user.

Properties:
$TYPEDFIELDS
"""
struct SoilStratigraphy{NF, Horizons <: NamedTuple} <: AbstractStratigraphy{NF}
    "Named tuple of soil horizons ordered from top to bottom"
    horizons::Horizons
end

function SoilStratigraphy(::Type{NF}; horizons...) where {NF}
    @assert length(horizons) > 0 "At least one soil horizon must be specified; for simple configurations, consider using HomogeneousSoilStratigraphy."
    horizon_nt = (; horizons...)
    return SoilStratigraphy{NF, typeof(horizon_nt)}(horizon_nt)
end

# Convenience constructors

"""
    $TYPEDSIGNATURES

Convenience constructor that creates a `SoilStratigraphy` with a single `ConstantSoilHorizon`.
"""
function HomogeneousSoilStratigraphy(
        ::Type{NF};
        texture::SoilTexture = SoilTexture(NF),
        porosity::AbstractSoilPorosity = ConstantSoilPorosity(NF)
    ) where {NF}
    return SoilStratigraphy(NF, horizon = ConstantSoilHorizon(NF; texture, porosity))
end

"""
    $TYPEDSIGNATURES

Convenience constructor that creates a three-layer `SoilStratigraphy` with spatially varying
organic (O), surface (A), and bedrock (R) horizons.
"""
function OARSoilStratigraphy(
        ::Type{NF};
        organic::AbstractSoilHorizon = PrescribedSoilHorizon(NF),
        surface::AbstractSoilHorizon = PrescribedSoilHorizon(NF),
        bedrock::AbstractSoilHorizon = PrescribedSoilHorizon(NF)
    ) where {NF}
    return SoilStratigraphy(NF; organic, surface, bedrock)
end

"""
    $TYPEDSIGNATURES

Convenience constructor that creates a three-layer `SoilStratigraphy` with spatially varying
organic (O), surface (A), subsoil (B), substratum (C), and bedrock (R) horizons.
"""
function OABCRSoilStratigraphy(
        ::Type{NF};
        organic::AbstractSoilHorizon = PrescribedSoilHorizon(NF),
        surface::AbstractSoilHorizon = PrescribedSoilHorizon(NF),
        subsoil::AbstractSoilHorizon = PrescribedSoilHorizon(NF),
        substratum::AbstractSoilHorizon = PrescribedSoilHorizon(NF),
        bedrock::AbstractSoilHorizon = PrescribedSoilHorizon(NF)
    ) where {NF}
    return SoilStratigraphy(NF; organic, surface, subsoil, substratum, bedrock)
end

Base.length(strat::SoilStratigraphy) = length(strat.horizons)
Base.iterate(strat::SoilStratigraphy) = Base.iterate(strat.horizons)
Base.iterate(strat::SoilStratigraphy, iter) = Base.iterate(strat.horizons, iter)

# Variables

function variables(strat::SoilStratigraphy)
    horizon_vars = map(variables, strat.horizons)
    return namespaces(horizon_vars)
end

# Kernel functions

"""
    $TYPEDSIGNATURES

Determine the index of the soil horizon in `strat` which contains the grid cell at index `i, j, k`.
"""
@inline function soil_horizon_index(i, j, k, grid, fields, strat::SoilStratigraphy)
    fgrid = get_field_grid(grid)
    # get midpoint of current node at k
    zₖ = znode(i, j, k, fgrid, Center(), Center(), Center())
    # initially set z to uppermost layer boundary
    z = znode(i, j, fgrid.Nz + 1, fgrid, Center(), Center(), Face())
    iₖ = 1
    for (l, name) in enumerate(keys(strat.horizons))
        horizon = strat.horizons[name]
        hfields = horizon_fields(fields, strat, l)
        iₖ = ifelse(zₖ <= z, l, iₖ)
        Δzₗ = soil_thickness(i, j, grid, hfields, horizon)
        z -= Δzₗ
    end
    return iₖ
end

"""
    $SIGNATURES

Retrieve the namespaced fields belonging to the `l`-th horizon of `strat` from `fields`.
"""
@inline horizon_fields(fields, strat::SoilStratigraphy, l::Integer) = fields.namespaces[keys(strat.horizons)[l]]

@inline function soil_horizon(i, j, k, grid, fields, strat::SoilStratigraphy)
    l = soil_horizon_index(i, j, k, grid, fields, strat)
    return strat.horizons[l]
end

@inline function soil_texture(i, j, k, grid, fields, strat::SoilStratigraphy)
    l = soil_horizon_index(i, j, k, grid, fields, strat)
    horizon = strat.horizons[l]
    texture = soil_texture(i, j, grid, horizon_fields(fields, strat, l), horizon)
    return texture
end

"""
    $TYPEDSIGNATURES

Compute the organic fraction of solid material in the soil volume at index `i, j, k`.
"""
@inline function organic_fraction(
        i, j, k, grid, fields,
        strat::SoilStratigraphy,
        bgc::AbstractSoilBiogeochemistry
    )
    ρ_soc = density_soc(i, j, k, grid, fields, bgc)
    ρ_org = density_pure_soc(bgc)
    l = soil_horizon_index(i, j, k, grid, fields, strat)
    horizon = strat.horizons[l]
    texture = soil_texture(i, j, grid, horizon_fields(fields, strat, l), horizon)
    por = organic_porosity(i, j, k, grid, fields, horizon.porosity, texture)
    organic = ρ_soc / ((1 - por) * ρ_org)
    return organic
end

"""
    $TYPEDSIGNATURES

Compute the porosity of the soil volume at the given indices.
"""
@propagate_inbounds function porosity(
        i, j, k, grid, fields,
        strat::SoilStratigraphy,
        bgc::AbstractSoilBiogeochemistry
    )
    organic = organic_fraction(i, j, k, grid, fields, strat, bgc)
    l = soil_horizon_index(i, j, k, grid, fields, strat)
    horizon = strat.horizons[l]
    texture = soil_texture(i, j, grid, horizon_fields(fields, strat, l), horizon)
    por_m = mineral_porosity(i, j, k, grid, fields, horizon.porosity, texture)
    por_o = organic_porosity(i, j, k, grid, fields, horizon.porosity, texture)
    return (1 - organic) * por_m + organic * por_o
end

"""
    $SIGNATURES

Construct a `SoilVolume` object summarizing the material composition of the soil volume
at the given indices `i, j, k` on `grid`.
"""
@propagate_inbounds function soil_volume(
        i, j, k, grid, fields,
        strat::AbstractStratigraphy,
        hydrology::AbstractSoilHydrology,
        bgc::AbstractSoilBiogeochemistry
    )
    # get current saturation and liquid fraction
    sat = saturation_water_ice(i, j, k, grid, fields, hydrology)
    liq = liquid_water_fraction(i, j, k, grid, fields, hydrology)
    # compute porosity and solid fractions
    por = porosity(i, j, k, grid, fields, strat, bgc)
    solid = soil_matrix(i, j, k, grid, fields, strat, bgc)
    return SoilVolume(por, sat, liq, solid)
end

"""
    $TYPEDSIGNATURES

Compute and return the soil solid matrix at index `i, j, k` on `grid`. The default
implementation assumes a simple `MineralOrganic` parameterization of the solid material.
"""
@propagate_inbounds function soil_matrix(
        i, j, k, grid, fields,
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry
    )
    organic = organic_fraction(i, j, k, grid, fields, strat, bgc)
    texture = soil_texture(i, j, k, grid, fields, strat)
    return MineralOrganic(texture, organic)
end
