struct SurfaceHydrology{
    NF,
    CanopyHydrology<:AbstractCanopyHydrology,
    Evapotranspiration<:AbstractEvapotranspiration,
    SurfaceRunoff<:AbstractSurfaceRunoff,
} <: AbstractSurfaceHydrology
    "Canopy hydrology scheme"
    canopy_hydrology::CanopyHydrology

    "Canopy evapotranspiration scheme"
    evapotranpsiration::Evapotranspiration
    
    "Surface runoff scheme"
    surface_runoff::SurfaceRunoff
end

function SurfaceHydrology(
    ::Type{NF};
    canopy_hydrology = PALADYNCanopyHydrology(NF),
    canopy_ET = PALADYNCanopyEvapotranspiration(NF),
    surface_runoff = DirectSurfaceRunoff(NF)
) where {NF}
    return SurfaceHydrology{NF, typeof(canopy_hydrology), typeof(canopy_ET), typeof(surface_runoff)}(canopy_hydrology, canopy_ET, surface_runoff)
end

variables(hydrology::SurfaceHydrology) = tuplejoin(
    variables(hydrology.canopy_hydrology),
    variables(hydrology.evapotranpsiration),
    variables(hydrology.surface_runoff)
)

function compute_auxiliary!(state, model, hydrology::SurfaceHydrology)
    compute_auxiliary!(state, model, hydrology.canopy_hydrology)
    compute_auxiliary!(state, model, hydrology.evapotranpsiration)
    compute_auxiliary!(state, model, hydrology.surface_runoff)
end

function compute_tendencies!(state, model, hydrology::SurfaceHydrology)
    compute_auxiliary!(state, model, hydrology.canopy_hydrology)
    compute_auxiliary!(state, model, hydrology.evapotranpsiration)
    compute_auxiliary!(state, model, hydrology.surface_runoff)
end
