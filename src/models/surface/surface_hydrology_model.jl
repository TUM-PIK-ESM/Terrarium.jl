"""
    $TYPEDEF

Model for surface hydrology processes.

Properties:
$TYPEDFIELDS
"""
@kwdef struct SurfaceHydrologyModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Atmosphere<:AbstractAtmosphere,
    CanopyHydrology<:AbstractCanopyHydrology,
    CanopyET<:AbstractEvapotranspiration,
    SurfaceRunoff<:AbstractSurfaceRunoff,
    Constants<:PhysicalConstants{NF},
    Initializer<:AbstractInitializer,
} <: AbstractSurfaceHydrologyModel{NF, GridType}
    "Spatial grid type"
    grid::GridType

    "Atmospheric input configuration"
    atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Canopy hydrology scheme"
    canopy_hydrology::CanopyHydrology = PALADYNCanopyHydrology(eltype(grid))

    "Canopy evapotranspiration scheme"
    evapotranpsiration::CanopyET = PALADYNCanopyEvapotranspiration(eltype(grid))
    
    "Surface runoff scheme"
    surface_runoff::SurfaceRunoff = PALADYNSurfaceRunoff(eltype(grid))

    "Physical constants"
    constants::Constants = PhysicalConstants(eltype(grid))

    "State variable initializer"
    initializer::Initializer = DefaultInitializer()
end

# Getter methods
get_grid(model::SurfaceHydrologyModel) = model.grid

get_atmosphere(model::SurfaceHydrologyModel) = model.atmosphere

get_surface_hydrology(model::SurfaceHydrologyModel) = SurfaceHydrology(
    model.canopy_hydrology,
    model.evapotranpsiration,
    model.surface_runoff
)

get_constants(model::SurfaceHydrologyModel) = model.constants

# Model interface methods
variables(model::SurfaceHydrologyModel) = tuplejoin(
    variables(model.atmosphere),
    variables(model.canopy_hydrology),
    variables(model.evapotranpsiration),
    variables(model.surface_runoff),
)

processes(model::SurfaceHydrologyModel) = (
    model.atmosphere,
    model.canopy_hydrology,
    model.evapotranpsiration,
    model.surface_runoff
)

function compute_auxiliary!(state, model::SurfaceHydrologyModel)
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.canopy_hydrology)
    compute_auxiliary!(state, model, model.evapotranpsiration)
    compute_auxiliary!(state, model, model.surface_runoff)
end

function compute_tendencies!(state, model::SurfaceHydrologyModel)
    compute_tendencies!(state, model, model.atmosphere)
    compute_tendencies!(state, model, model.canopy_hydrology)
    compute_tendencies!(state, model, model.evapotranpsiration)
    compute_tendencies!(state, model, model.surface_runoff)
end

