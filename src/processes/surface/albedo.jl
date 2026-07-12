"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
@kwdef struct PrescribedAlbedo{NF} <: AbstractAlbedo{NF} end

PrescribedAlbedo(::Type{NF}) where {NF} = PrescribedAlbedo{NF}()

variables(::PrescribedAlbedo) = (
    input(:albedo, XY(), bounds = UnitInterval, desc = "Surface albedo, i.e. ratio of outgoing to incoming shortwave radiation [-]"),
    input(:emissivity, XY(), bounds = UnitInterval, desc = "Surface emissivity, i.e. efficiency of longwave emission [-]"),
)

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
@parameterized @kwdef struct ConstantAlbedo{NF} <: AbstractAlbedo{NF}
    "Surface albedo, i.e. ratio of outgoing to incoming shortwave radiation"
    @param albedo::NF = 0.3 (bounds = UnitInterval,)

    "Surface emissivity, i.e. fraction of thermal radiation emitted from the surface"
    @param emissivity::NF = 0.97 (bounds = UnitInterval,)
end

ConstantAlbedo(::Type{NF}; kwargs...) where {NF} = ConstantAlbedo{NF}(; kwargs...)

"""
    $TYPEDSIGNATURES

Return the surface albedo at grid point `i, j`.
"""
@inline albedo(i, j, grid, fields, albedo::ConstantAlbedo) = albedo.albedo

"""
    $TYPEDSIGNATURES

Return the surface emissivity at grid point `i, j`.
"""
@inline emissivity(i, j, grid, fields, albedo::ConstantAlbedo) = albedo.emissivity

"""
    $TYPEDEF

Diagnosed surface albedo and emissivity. The values are computed in [`compute_auxiliary!`](@ref) as a
snow-cover-weighted blend of a snow-free background and snow, `α = (1 − f_snow)·α_bg + f_snow·α_snow`
(and likewise for emissivity), where `f_snow` is the snow-covered area fraction of the optional snow
component passed to `compute_auxiliary!`. Without a snow component (`snow === nothing`), `f_snow = 0` and
the background values are used.

Properties:
$TYPEDFIELDS
"""
@parameterized @kwdef struct DiagnosticAlbedo{NF} <: AbstractAlbedo{NF}
    "Snow-free (background) albedo"
    @param background_albedo::NF = 0.25 (bounds = UnitInterval,)

    "Snow-free (background) emissivity"
    @param background_emissivity::NF = 0.97 (bounds = UnitInterval,)

    "Snow albedo"
    @param snow_albedo::NF = 0.8 (bounds = UnitInterval,)

    "Snow emissivity"
    @param snow_emissivity::NF = 0.99 (bounds = UnitInterval,)
end

DiagnosticAlbedo(::Type{NF}; kwargs...) where {NF} = DiagnosticAlbedo{NF}(; kwargs...)

variables(::DiagnosticAlbedo) = (
    auxiliary(:albedo, XY(), bounds = UnitInterval, desc = "Diagnosed surface albedo [-]"),
    auxiliary(:emissivity, XY(), bounds = UnitInterval, desc = "Diagnosed surface emissivity [-]"),
)

"""
    $TYPEDSIGNATURES

Diagnose the blended surface albedo and emissivity, weighting the snow-free background and snow values
by the snow-covered area fraction of the optional `snow` component.
"""
function compute_auxiliary!(
        state, grid,
        alb::DiagnosticAlbedo,
        snow = nothing,
        args...
    )
    out = auxiliary_fields(state, alb)
    # include the snow fields (for `snow_cover_fraction`) only when a snow component is present
    fields = isnothing(snow) ? get_fields(state, alb; except = out) : get_fields(state, alb, snow; except = out)
    launch!(grid, XY, compute_albedo_kernel!, out, fields, alb, snow)
    return nothing
end

# Kernel functions

""" $TYPEDSIGNATURES """
@propagate_inbounds function compute_albedo!(out, i, j, grid, fields, alb::DiagnosticAlbedo, snow)
    f = snow_cover_fraction(i, j, grid, fields, snow)
    out.albedo[i, j, 1] = (one(f) - f) * alb.background_albedo + f * alb.snow_albedo
    out.emissivity[i, j, 1] = (one(f) - f) * alb.background_emissivity + f * alb.snow_emissivity
    return nothing
end

# Kernels

@kernel inbounds = true function compute_albedo_kernel!(out, grid, fields, alb::DiagnosticAlbedo, snow)
    i, j = @index(Global, NTuple)
    compute_albedo!(out, i, j, grid, fields, alb, snow)
end
