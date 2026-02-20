"""
    $TYPEDEF

Simple surface runoff scheme that computes runoff as

```math
R = P + D - I
```
where `P` is precipitation reaching the ground, `D` is drainage from accumualted excess
water at the surface, and `I` is infiltration into the soil.

Properties:
$FIELDS
"""
@kwdef struct DirectSurfaceRunoff{NF} <: AbstractSurfaceRunoff{NF}
    "Surface water removal timescale"
    τ_r::NF = 3600.0
end

DirectSurfaceRunoff(::Type{NF}; kwargs...) where {NF} = DirectSurfaceRunoff{NF}(; kwargs...)

"""
    $TYPEDSIGNATURES

Compute surface drainage flux from the current `surface_excess_water` resevoir state.
"""
@inline function compute_surface_drainage(runoff::DirectSurfaceRunoff{NF}, surface_excess_water) where {NF}
    let S = max(surface_excess_water, zero(NF)),
            τ = runoff.τ_r
        ∂S∂t = S / τ
        return ∂S∂t
    end
end

"""
    $TYPEDSIGNATURES

Compute infiltration from the given `influx` (water available for infiltration), saturation of the uppermost
soil layer `sat_top`, and the maximum allowed infiltration `max_infil`.
"""
@inline function compute_infiltration(runoff::DirectSurfaceRunoff{NF}, influx, sat_top, max_infil) where {NF}
    let is_unsaturated = sat_top < one(NF)
        # Infiltration is min of
        infil = min(influx, max_infil) * is_unsaturated
        return infil
    end
end

"""
    $TYPEDSIGNATURES

Compute surface runoff as `precipitation + surface_drainage - infiltration`.
"""
@inline function compute_surface_runoff(runoff::DirectSurfaceRunoff, precip_ground, surface_drainage, infil)
    let P = precip_ground,
            ∂S∂t = surface_drainage,
            I = infil
        # Compute runoff as residual of precipitation + drainage - infiltration
        surface_runoff = P + ∂S∂t - I
        return surface_runoff
    end
end

# Process methods

variables(::DirectSurfaceRunoff) = (
    auxiliary(:surface_runoff, XY(), units = u"m/s", desc = "Total surface runoff"),
    auxiliary(:infiltration, XY(), units = u"m/s", desc = "Infiltration flux"),
)

function compute_auxiliary!(
        state, grid,
        runoff::DirectSurfaceRunoff,
        canopy_interception::AbstractCanopyInterception,
        soil::AbstractSoil,
        args...
    )
    soil_hydrology = get_hydrology(soil)
    out = auxiliary_fields(state, runoff)
    fields = get_fields(state, runoff, canopy_interception, soil_hydrology; except = out)
    launch!(grid, XY, compute_auxiliary_kernel!, out, fields, runoff, canopy_interception, soil_hydrology)
    return nothing
end

# Kernel function

@propagate_inbounds function compute_surface_runoff!(
        out, i, j, grid, fields,
        runoff::DirectSurfaceRunoff{NF},
        canopy_interception::AbstractCanopyInterception,
        soil_hydrology::AbstractSoilHydrology
    ) where {NF}
    fgrid = get_field_grid(grid)

    # Get inputs
    precip_ground = ground_precipitation(i, j, grid, fields, canopy_interception)
    excess_water = surface_excess_water(i, j, grid, fields, soil_hydrology)
    k_unsat = hydraulic_conductivity(i, j, fgrid.Nz, grid, fields, soil_hydrology)
    sat_top = saturation_water_ice(i, j, fgrid.Nz, grid, fields, soil_hydrology)

    # Case 1: Excess water present at the surface -> precipitation adds to excess water
    # and we set the infiltration rate to the min of hydraulic conductivity and surface_excess_water
    if excess_water > zero(NF)
        # Compute rate of excess water removal (surface drainage)
        surface_drainage = compute_surface_drainage(runoff, excess_water)
        # Calculate infiltration
        infil = out.infiltration[i, j, 1] = compute_infiltration(runoff, surface_drainage, sat_top, k_unsat)
        # Case 2: No excess water -> rainfall is routed directly to infiltration
    else
        surface_drainage = zero(NF)
        infil = out.infiltration[i, j, 1] = compute_infiltration(runoff, precip_ground, sat_top, k_unsat)
    end

    # Compute surface runoff
    out.surface_runoff[i, j, 1] = compute_surface_runoff(runoff, precip_ground, surface_drainage, infil)
    return out
end

# Kernels

@kernel inbounds = true function compute_auxiliary_kernel!(out, grid, fields, runoff::AbstractSurfaceRunoff, args...)
    i, j = @index(Global, NTuple)
    compute_surface_runoff!(out, i, j, grid, fields, runoff, args...)
end
