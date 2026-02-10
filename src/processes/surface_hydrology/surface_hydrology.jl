"""
    $TYPEDEF

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
        canopy_interception = PALADYNCanopyInterception(NF),
        canopy_ET = PALADYNCanopyEvapotranspiration(NF),
        surface_runoff = DirectSurfaceRunoff(NF)
    ) where {NF}
    return SurfaceHydrology{NF, typeof(canopy_interception), typeof(canopy_ET), typeof(surface_runoff)}(canopy_interception, canopy_ET, surface_runoff)
end

function compute_auxiliary!(
        state, grid,
        hydrology::SurfaceHydrology,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
        soil::AbstractSoil
    )
    compute_auxiliary!(state, grid, hydrology.canopy_interception, atmos, constants)
    compute_auxiliary!(state, grid, hydrology.evapotranspiration, hydrology.canopy_interception, atmos, constants, soil)
    compute_auxiliary!(state, grid, hydrology.surface_runoff, hydrology.canopy_interception, soil)
    return nothing
end

function compute_tendencies!(
        state, grid,
        hydrology::SurfaceHydrology,
        args...,
    )
    # Compute tendencies for canopy interception
    compute_tendencies!(state, grid, hydrology.canopy_interception, hydrology.evapotranspiration)
    return nothing
end
