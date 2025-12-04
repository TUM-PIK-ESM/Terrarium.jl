"""
    $TYPEDEF

Model for surface hydrology processes.

Properties:
$TYPEDFIELDS
"""
@kwdef struct SurfaceHydrologyModel{
    NF,
    CanopyHydrology<:AbstractCanopyHydrology,
    Runoff<:AbstractRunoff,
    Snow<:AbstractSnow,
    GridType<:AbstractLandGrid{NF},
    Atmosphere<:AbstractAtmosphere,
    Initializer<:AbstractInitializer,
} <: AbstractSurfaceHydrologyModel{NF, GridType}
    "Canopy hydrology scheme"
    canopy_hydrology::CanopyHydrology = CanopyHydrology()
    
    "Runoff scheme"
    runoff::Runoff = Runoff()

    "Snow scheme"
    snow::Snow = Snow()

    "Spatial grid type"
    grid::GridType

    "Atmospheric inputs"
    atmosphere::Atmosphere

    "Physical constants"
    constants::PhysicalConstants

    "State variable initializer"
    initializer::Initializer
end

# SurfaceHydrologyModel getter methods
get_atmosphere(model::SurfaceHydrologyModel) = model.atmosphere

get_canopy_hydrology(model::SurfaceHydrologyModel) = model.canopy_hydrology

get_runoff(model::SurfaceHydrologyModel) = model.runoff

get_snow(model::SurfaceHydrologyModel) = model.snow

get_constants(model::SurfaceHydrologyModel) = model.constants

# Model interface methods
variables(model::SurfaceHydrologyModel) = tuplejoin(
    variables(model.atmosphere)...,
    variables(model.canopy_hydrology)...,
    variables(model.runoff)...,
    variables(model.snow)...,
)

get_processes(model::SurfaceHydrologyModel) = (
    model.atmosphere,
    model.canopy_hydrology,
    model.runoff,
    model.snow,
)

function compute_auxiliary!(state, model::SurfaceHydrologyModel)
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.canopy_hydrology)
    compute_auxiliary!(state, model, model.runoff)
    compute_auxiliary!(state, model, model.snow)
end

function compute_tendencies!(state, ::SurfaceHydrologyModel)
    compute_tendencies!(state, model, model.atmosphere)
    compute_tendencies!(state, model, model.canopy_hydrology)
    compute_tendencies!(state, model, model.runoff)
    compute_tendencies!(state, model, model.snow)
end

