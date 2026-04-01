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
        canopy_interception::CI = PALADYNCanopyInterception(NF),
        evapotranspiration::ET = PALADYNCanopyEvapotranspiration(NF),
        surface_runoff::SR = DirectSurfaceRunoff(NF)
    ) where {NF, CI, ET, SR}
    return SurfaceHydrology{NF, CI, ET, SR}(canopy_interception, evapotranspiration, surface_runoff)
end

function compute_auxiliary!(
        state, grid,
        hydrology::SurfaceHydrology,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
        soil::Optional{AbstractSoil} = nothing,
        vegetation::Optional{AbstractVegetation} = nothing,
        args...
    )
    compute_auxiliary!(state, grid, hydrology.canopy_interception, atmos, constants)
    compute_auxiliary!(state, grid, hydrology.evapotranspiration, atmos, constants, soil, vegetation, hydrology.canopy_interception)
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
