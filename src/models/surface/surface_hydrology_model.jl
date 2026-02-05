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
    CanopyHydrology<:AbstractCanopyInterception,
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
    canopy_interception::CanopyHydrology = PALADYNCanopyInterception(eltype(grid))

    "Canopy evapotranspiration scheme"
    evapotranpsiration::CanopyET = PALADYNCanopyEvapotranspiration(eltype(grid))
    
    "Surface runoff scheme"
    surface_runoff::SurfaceRunoff = DirectSurfaceRunoff(eltype(grid))

    "Physical constants"
    constants::Constants = PhysicalConstants(eltype(grid))

    "State variable initializer"
    initializer::Initializer = DefaultInitializer()
end

# Model interface methods

function compute_auxiliary!(state, model::SurfaceHydrologyModel)
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.canopy_interception)
    compute_auxiliary!(state, model, model.evapotranpsiration)
    compute_auxiliary!(state, model, model.surface_runoff)
end

function compute_tendencies!(state, model::SurfaceHydrologyModel)
    compute_tendencies!(state, model, model.atmosphere)
    compute_tendencies!(state, model, model.canopy_interception)
    compute_tendencies!(state, model, model.evapotranpsiration)
    compute_tendencies!(state, model, model.surface_runoff)
end

