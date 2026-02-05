"""
    $TYPEDEF

Properties:
$FIELDS
"""
struct SurfaceHydrology{
    NF,
    CanopyHydrology<:AbstractCanopyInterception{NF},
    Evapotranspiration<:AbstractEvapotranspiration{NF},
    SurfaceRunoff<:AbstractSurfaceRunoff{NF},
} <: AbstractSurfaceHydrology{NF}
    "Canopy hydrology scheme"
    canopy_hydrology::CanopyHydrology

    "Canopy evapotranspiration scheme"
    evapotranpsiration::Evapotranspiration
    
    "Surface runoff scheme"
    surface_runoff::SurfaceRunoff
end

function SurfaceHydrology(
    ::Type{NF};
    canopy_hydrology = PALADYNCanopyInterception(NF),
    canopy_ET = PALADYNCanopyEvapotranspiration(NF),
    surface_runoff = DirectSurfaceRunoff(NF)
) where {NF}
    return SurfaceHydrology{NF, typeof(canopy_hydrology), typeof(canopy_ET), typeof(surface_runoff)}(canopy_hydrology, canopy_ET, surface_runoff)
end

function compute_auxiliary!(
    state, grid,
    hydrology::SurfaceHydrology,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    compute_auxiliary!(state, grid, hydrology.canopy_hydrology, atmos, constants)
    compute_auxiliary!(state, grid, hydrology.evapotranpsiration)
    compute_auxiliary!(state, grid, hydrology.surface_runoff)
end

function compute_tendencies!(
    state, grid,
    hydrology::SurfaceHydrology,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    compute_auxiliary!(state, grid, hydrology.canopy_hydrology)
    compute_auxiliary!(state, grid, hydrology.evapotranpsiration)
    compute_auxiliary!(state, grid, hydrology.surface_runoff)
end
