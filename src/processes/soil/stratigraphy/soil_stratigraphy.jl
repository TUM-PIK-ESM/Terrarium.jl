"""
    $TYPEDEF

Represents a soil stratigraphy as a stack of named [soil horizons](https://en.wikipedia.org/wiki/Soil_horizon).
Each soil horizon is assumed to have internally homogeneous soil properties. The number of horizons and their
respective names are defined by the user.

Properties:
$TYPEDFIELDS
"""
struct SoilStratigraphy{NF, N, Horizons <: Tuple{Vararg{AbstractSoilHorizon{NF}, N}}} <: AbstractStratigraphy{NF}
    "Named tuple of soil horizons ordered from top to bottom"
    horizons::Horizons
end

function SoilStratigraphy(::Type{NF}, horizons::AbstractSoilHorizon...) where {NF}
    @assert length(horizons) > 0 "At least one soil horizon must be specified; for simple configurations, consider using HomogeneousSoilStratigraphy."
    return SoilStratigraphy(horizons)
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
    return SoilStratigraphy(NF, ConstantSoilHorizon(NF, :soil; texture, porosity))
end

Base.length(strat::SoilStratigraphy) = length(strat.horizons)
Base.iterate(strat::SoilStratigraphy) = Base.iterate(strat.horizons)
Base.iterate(strat::SoilStratigraphy, iter) = Base.iterate(strat.horizons, iter)

# Variables

function variables(strat::SoilStratigraphy)
    return map(strat.horizons) do horizon
        namespace(nameof(horizon), variables(horizon))
    end
end

# Process methods

compute_auxiliary!(state, grid, ::SoilStratigraphy, args...) = nothing

compute_tendencies!(state, grid, ::SoilStratigraphy, args...) = nothing

# Kernel functions

"""
    $TYPEDSIGNATURES

Retrieve the soil horizon for the soil volume at index `i, j, k`. The last
soil horizon in `strat` is assumed to extend to the bottom of the vertical
column regardless of its associated `thickness`. Note that, since `SoilStratigraphy`
uses namespaces for the state variables of each horizon, any methods defined on
`AbstractSoilHorizon` types should be passed the namespace, e.g:
```julia
horizon = soil_horizon(i, j, k, grid, fields, strat)
texture = soil_texture(i, j, grid, getproperty(fields, nameof(horizon)), horizon)
```
"""
@inline function soil_horizon(i, j, k, grid, fields, strat::SoilStratigraphy{NF}) where {NF}
    fgrid = get_field_grid(grid)
    # get midpoint of current node at k
    zₖ = znode(i, j, k, fgrid, Center(), Center(), Center())
    # initially set z to uppermost layer boundary
    z = znode(i, j, fgrid.Nz + 1, fgrid, Center(), Center(), Face())
    # initialize result to first horizon (will be updated in loop)
    horizon = first(strat.horizons)
    for next in strat.horizons
        # update result when zₖ is within this horizon's depth range
        horizon = ifelse(zₖ <= z, next, horizon)
        Δzₗ = soil_thickness(i, j, grid, getproperty(fields, nameof(next)), next)
        z -= Δzₗ
    end
    return horizon
end

"""
    $TYPEDSIGNATURES

Convenience method that invokes `func(i, j, grid, horizon_fields, horizon, args...; kwargs...)`
where `horizon` is the soil horizon returned by [`soil_horizon`](@ref) and `horizon_fields`
is `getproperty(fields, nameof(horizon))`.
"""
@inline function with_soil_horizon(
        func, i, j, k, grid, fields,
        strat::SoilStratigraphy,
        args...;
        kwargs...
    )
    horizon = soil_horizon(i, j, k, grid, fields, strat)
    return func(i, j, grid, getproperty(fields, nameof(horizon)), horizon, args...; kwargs...)
end

"""
    $TYPEDSIGNATURES

Retrieve the soil texture of the soil volume at index `i, j, k` in the given stratigraphy `strat`.
"""
@inline function soil_texture(i, j, k, grid, fields, strat::SoilStratigraphy{NF}) where {NF}
    texture = with_soil_horizon(soil_texture, i, j, k, grid, fields, strat)
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
    por = with_soil_horizon(organic_porosity, i, j, k, grid, fields, strat)
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
    por_m = with_soil_horizon(mineral_porosity, i, j, k, grid, fields, strat)
    por_o = with_soil_horizon(organic_porosity, i, j, k, grid, fields, strat)
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
