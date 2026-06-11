"""
    $TYPEDEF

Represents a soil stratigraphy.

Properties:
$TYPEDFIELDS
"""
struct SoilStratigraphy{NF, Horizons <: NamedTuple} <: AbstractStratigraphy{NF}
    horizons::Horizons
end

function SoilStratigraphy(::Type{NF}; horizons...) where {NF}
    horizon_nt = (; horizons...)
    return SoilStratigraphy{NF, typeof(horizon_nt)}(horizon_nt)
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

@inline function soil_horizon(i, j, k, grid, fields, strat::SoilStratigraphy{NF}) where {NF}
    fgrid = get_field_grid(grid)
    # get midpoint of current node at k
    zₖ = znode(i, j, k, fgrid, Center(), Center(), Center())
    # initially set z to uppermost layer boundary
    z = znode(i, j, fgrid.Nz + 1, fgrid, Center(), Center(), Face())
    iₖ = 1
    for (l, name) in enumerate(keys(strat.horizons))
        horizon = strat.horizons[name]
        horizon_fields = fields.namespaces[name]
        iₖ = ifelse(zₖ <= z, l, iₖ)
        Δzₗ = soil_thickness(i, j, grid, horizon_fields, horizon)
        z -= Δzₗ
    end
    return strat.horizons[iₖ]
end

@inline function soil_texture(i, j, k, grid, fields, strat::SoilStratigraphy{NF}) where {NF}
    horizon = soil_horizon(i, j, k, grid, fields, strat)
    texture = soil_texture(i, j, grid, fields, horizon)
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
    horizon = soil_horizon(i, j, k, grid, fields, strat)
    por = organic_porosity(i, j, k, grid, fields, horizon.porosity, horizon.texture)
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
    horizon = soil_horizon(i, j, k, grid, fields, strat)
    por_m = mineral_porosity(i, j, k, grid, fields, horizon.porosity, horizon.texture)
    por_o = organic_porosity(i, j, k, grid, fields, horizon.porosity, horizon.texture)
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
