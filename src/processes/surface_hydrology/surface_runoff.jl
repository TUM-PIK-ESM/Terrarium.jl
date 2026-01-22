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
@kwdef struct DirectSurfaceRunoff{NF} <: AbstractSurfaceRunoff
    "Surface water removal timescale"
    τ_r::NF = 3600.0
end

DirectSurfaceRunoff(::Type{NF}; kwargs...) where {NF} = DirectSurfaceRunoff{NF}(; kwargs...)

variables(::DirectSurfaceRunoff) = (
    auxiliary(:surface_runoff, XY(), units=u"m/s", desc="Total surface runoff"),
    auxiliary(:infiltration, XY(), units=u"m/s", desc="Infiltration flux"),
)

function compute_auxiliary!(state, model, runoff::DirectSurfaceRunoff)
    grid = get_grid(model)
    soil_hydrology = get_soil_hydrology(model)
    surface_hydrology = get_surface_hydrology(model)
    launch!(state, grid, :xy, compute_auxiliary_kernel!, runoff, surface_hydrology.canopy_hydrology, soil_hydrology)
end

# Kernels

@kernel function compute_auxiliary_kernel!(
    state, grid,
    runoff::DirectSurfaceRunoff{NF},
    canopy_hydrology::AbstractCanopyHydrology,
    soil_hydrology::AbstractSoilHydrology
) where {NF}
    i, j = @index(Global, NTuple)
    fgrid = get_field_grid(grid)

    @inbounds let precip_ground = ground_precipitation(i, j, state, grid, canopy_hydrology),
                  surface_excess_water = surface_excess_water(i, j, state, grid, soil_hydrology),
                  k_unsat = hydraulic_conductivity(i, j, fgrid.Nz, state, grid, soil_hydrology),
                  sat_top = saturation_water_ice(i, j, fgrid.Nz, state, grid, soil_hydrology);
        # Case 1: Excess water present at the surface -> precipitation adds to excess water
        # and we set the infiltration rate to the min of hydraulic conductivity and surface_excess_water
        if surface_excess_water > zero(NF)
            # Compute rate of excess water removal (surface drainage)
            surface_drainage = compute_surface_drainage(runoff, surface_excess_water)
            # Calculate infiltration
            infil = state.infiltration[i, j, 1] = compute_infiltration(runoff, sat_top, k_unsat, surface_drainage)
        # Case 2: No excess water -> rainfall is routed directly to infiltration
        else
            surface_drainage = zero(NF)
            infil = state.infiltration[i, j, 1] = min(k_unsat, precip_ground) * (sat_top < one(NF))
        end
        # Compute runoff as residual of precipitation + drainage - infiltration
        state.surface_runoff[i, j, 1] = compute_surface_runoff(runoff, precip_ground, surface_drainage, infil)
    end
end

# Kernel functions

@inline function compute_surface_drainage(runoff::DirectSurfaceRunoff, surface_excess_water)
    let S = surface_excess_water,
        τ = runoff.τ_r;
        ∂S∂t = S / τ
        return ∂S∂t
    end
end

@inline function compute_infiltration(runoff::DirectSurfaceRunoff{NF}, sat_top, k_unsat, surface_drainage) where {NF}
    let is_unsaturated = sat_top < one(NF),
        K = k_unsat,
        ∂S∂t = surface_drainage;
        infil = min(K, ∂S∂t) * is_unsaturated
        return infil
    end
end

@inline function compute_surface_runoff(runoff::DirectSurfaceRunoff, precip_ground, surface_drainage, infil)
    let P = precip_ground,
        ∂S∂t = surface_drainage,
        I = infil;
        surface_runoff = P + ∂S∂t - I
        return surface_runoff
    end
end
