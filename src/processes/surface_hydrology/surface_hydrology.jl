"""
    $TYPEDEF

Default representation of coupled surface hydrology processes including
canopy rain/snow interception, evapotranspiration, and surface runoff.

Properties:
$FIELDS
"""
struct SurfaceHydrology{
        NF,
        CanopyInterception <: AbstractCanopyInterception{NF},
        Evapotranspiration <: AbstractEvapotranspiration{NF},
        SurfaceRunoff <: AbstractSurfaceRunoff{NF},
    } <: AbstractSurfaceHydrology{NF}
    "Canopy hydrology scheme"
    canopy_interception::CanopyInterception

    "Canopy evapotranspiration scheme"
    evapotranspiration::Evapotranspiration

    "Surface runoff scheme"
    surface_runoff::SurfaceRunoff
end

function SurfaceHydrology(
        ::Type{NF};
        canopy_interception::CI = PALADYNCanopyInterception(NF),
        evapotranspiration::ET = PALADYNCanopyEvapotranspiration(NF),
        surface_runoff::SR = DirectSurfaceRunoff(NF)
    ) where {NF, CI, ET, SR}
    return SurfaceHydrology{NF, CI, ET, SR}(canopy_interception, evapotranspiration, surface_runoff)
end

""" $TYPEDSIGNATURES """
function compute_auxiliary!(
        state, grid,
        hydrology::SurfaceHydrology,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil} = nothing,
        vegetation::Optional{AbstractVegetation} = nothing,
        args...
    )
    compute_auxiliary!(state, grid, hydrology.canopy_interception, atmos)
    compute_auxiliary!(state, grid, hydrology.evapotranspiration, hydrology.canopy_interception, constants, atmos, soil, vegetation)
    compute_auxiliary!(state, grid, hydrology.surface_runoff, hydrology.canopy_interception, soil)
    return nothing
end

"""
    $TYPEDSIGNATURES

Recompute the evapotranspiration fluxes for `hydrology` from the current skin temperature by
re-running the ET auxiliary computation. Intended to be called after the surface energy balance
solve so that the partitioned ET fluxes (consumed by the soil-moisture forcing, canopy-water
tendency, and diagnostic output) are consistent with the converged skin temperature. The fluxes
the SEB itself consumes are evaluated lazily during the solve via `surface_humidity_flux` and do
not depend on this call.
"""
function compute_evapotranspiration!(
        state, grid,
        hydrology::SurfaceHydrology,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil} = nothing,
        vegetation::Optional{AbstractVegetation} = nothing,
        args...
    )
    compute_auxiliary!(state, grid, hydrology.evapotranspiration, hydrology.canopy_interception, constants, atmos, soil, vegetation)
    return nothing
end

""" $TYPEDSIGNATURES """
function compute_tendencies!(
        state, grid,
        hydrology::SurfaceHydrology,
        args...,
    )
    # Compute tendencies for canopy interception
    compute_tendencies!(state, grid, hydrology.canopy_interception, hydrology.evapotranspiration)
    return nothing
end
