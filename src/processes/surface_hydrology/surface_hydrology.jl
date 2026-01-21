@kwdef struct SurfaceHydrology{
    NF,
    Atmosphere<:AbstractAtmosphere,
    CanopyHydrology<:AbstractCanopyHydrology,
    Evapotranspiration<:AbstractEvapotranspiration,
    SurfaceRunoff<:AbstractSurfaceRunoff,
} <: AbstractSurfaceHydrology
    "Canopy hydrology scheme"
    canopy_hydrology::CanopyHydrology = PALADYNCanopyHydrology(eltype(grid))

    "Canopy evapotranspiration scheme"
    evapotranpsiration::Evapotranspiration = PALADYNCanopyEvapotranspiration(eltype(grid))
    
    "Surface runoff scheme"
    surface_runoff::SurfaceRunoff = PALADYNSurfaceRunoff(eltype(grid))
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
