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

Apply `f(i, j, grid, horizon_fields, horizon)` to the soil horizon containing the grid cell
at index `i, j, k`, where `horizon_fields` is the state variable namespace associated with
that horizon. The horizon stack is unrolled at compile time, so the call is type stable
(and thus GPU compatible) as long as `f` returns a value of the same type for every horizon.
"""
@inline function with_soil_horizon(f::F, i, j, k, grid, fields, strat::SoilStratigraphy) where {F}
    fgrid = get_field_grid(grid)
    # get midpoint of current node at k
    zₖ = znode(i, j, k, fgrid, Center(), Center(), Center())
    # initially set z to uppermost layer boundary
    z = znode(i, j, fgrid.Nz + 1, fgrid, Center(), Center(), Face())
    return _with_soil_horizon(f, i, j, zₖ, z, grid, fields, strat.horizons)
end

# Base case: the bottommost horizon contains all remaining nodes.
@inline function _with_soil_horizon(f::F, i, j, zₖ, z, grid, fields, horizons::NamedTuple{names, <:Tuple{Any}}) where {F, names}
    horizon_fields = getproperty(fields.namespaces, first(names))
    return f(i, j, grid, horizon_fields, first(horizons))
end

# Recursive case: apply f if zₖ lies above the lower boundary of the current horizon,
# otherwise recurse into the remaining horizons. Since the horizon names are encoded in
# the type, the recursion is fully unrolled by the compiler.
@inline function _with_soil_horizon(f::F, i, j, zₖ, z, grid, fields, horizons::NamedTuple{names}) where {F, names}
    horizon = first(horizons)
    horizon_fields = getproperty(fields.namespaces, first(names))
    Δz = soil_thickness(i, j, grid, horizon_fields, horizon)
    if zₖ > z - Δz
        return f(i, j, grid, horizon_fields, horizon)
    else
        return _with_soil_horizon(f, i, j, zₖ, z - Δz, grid, fields, Base.tail(horizons))
    end
end

"""
    $TYPEDSIGNATURES

Return the soil horizon containing the grid cell at index `i, j, k`. Note that when the
stratigraphy contains horizons of heterogeneous types, the return type depends on runtime
state, which makes this function unsuitable for use in GPU kernels; use
[`with_soil_horizon`](@ref) instead.
"""
@inline function soil_horizon(i, j, k, grid, fields, strat::SoilStratigraphy)
    return with_soil_horizon((i, j, grid, hfields, horizon) -> horizon, i, j, k, grid, fields, strat)
end

@inline function soil_texture(i, j, k, grid, fields, strat::SoilStratigraphy)
    return with_soil_horizon(soil_texture, i, j, k, grid, fields, strat)
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
    por = with_soil_horizon(i, j, k, grid, fields, strat) do i, j, grid, horizon_fields, horizon
        texture = soil_texture(i, j, grid, horizon_fields, horizon)
        organic_porosity(i, j, k, grid, horizon_fields, horizon.porosity, texture)
    end
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
    return with_soil_horizon(i, j, k, grid, fields, strat) do i, j, grid, horizon_fields, horizon
        texture = soil_texture(i, j, grid, horizon_fields, horizon)
        por_m = mineral_porosity(i, j, k, grid, horizon_fields, horizon.porosity, texture)
        por_o = organic_porosity(i, j, k, grid, horizon_fields, horizon.porosity, texture)
        (1 - organic) * por_m + organic * por_o
    end
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
