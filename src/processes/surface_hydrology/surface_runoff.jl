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
    auxiliary(:runoff_surface, XY(), units=u"m/s", desc="Surface runoff"),
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
                  surface_excess_water = state.surface_excess_water[i, j]
                  k_unsat = state.hydraulic_conductivity[i, j, fgrid.Nz],
                  sat_top = state.saturation_water_ice[i, j, 1];
        # Case 1: Excess water present at the surface -> precipitation adds to excess water
        # and we set the infiltration rate to the min of hydraulic conductivity and surface_excess_water
        if surface_excess_water > zero(NF)
            # Compute drainage of excess surface water
            surface_excess_water_drainage = surface_excess_water / runoff.τ_r
            # Calculate infiltration
            infil = state.infiltration[i, j, 1] = min(k_unsat, surface_excess_water_drainage) * (sat_top < one(NF))
            # Compute tendency for surface excess water
            state.tendencies.surface_excess_water[i, j, 1] = precip_ground - surface_excess_water_drainage
        # Case 2: No excess water -> rainfall is routed directly to infiltration
        else
            surface_excess_water_drainage = zero(NF)
            infil = state.infiltration[i, j, 1] = min(k_unsat, precip_ground) * (sat_top < one(NF))
        end
        # Compute runoff as residual of precipitation + drainage - infiltration
        state.runoff_surface[i, j, 1] = precip_ground + surface_excess_water_drainage - infil
    end
end
