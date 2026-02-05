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
    canopy_interception::CanopyHydrology

    "Canopy evapotranspiration scheme"
    evapotranpsiration::Evapotranspiration
    
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
    constants::PhysicalConstants
)
    compute_auxiliary!(state, grid, hydrology.canopy_interception, atmos, constants)
    compute_auxiliary!(state, grid, hydrology.evapotranpsiration)
    compute_auxiliary!(state, grid, hydrology.surface_runoff)
end

function compute_tendencies!(
    state, grid,
    hydrology::SurfaceHydrology,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    compute_tendencies!(state, grid, hydrology.canopy_interception)
    compute_tendencies!(state, grid, hydrology.evapotranpsiration)
    compute_tendencies!(state, grid, hydrology.surface_runoff)
end
