"""
    $TYPEDEF

Model for surface hydrology processes.

Properties:
$TYPEDFIELDS
"""
@kwdef struct SurfaceHydrologyModel{
        NF,
        GridType <: AbstractLandGrid{NF},
        Atmosphere <: AbstractAtmosphere,
        CanopyHydrology <: AbstractCanopyInterception,
        CanopyET <: AbstractEvapotranspiration,
        SurfaceRunoff <: AbstractSurfaceRunoff,
        Constants <: PhysicalConstants{NF},
        Initializer <: AbstractInitializer,
        Timesteppers <: NamedTuple,
    } <: AbstractSurfaceHydrologyModel{NF, GridType}
    "Spatial grid type"
    grid::GridType

    "Atmospheric input configuration"
    atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Canopy hydrology scheme"
    canopy_interception::CanopyHydrology = PALADYNCanopyInterception(eltype(grid))

    "Canopy evapotranspiration scheme"
    evapotranspiration::CanopyET = PALADYNCanopyEvapotranspiration(eltype(grid))

    "Surface runoff scheme"
    surface_runoff::SurfaceRunoff = DirectSurfaceRunoff(eltype(grid))

    "Physical constants"
    constants::Constants = PhysicalConstants(eltype(grid))

    "State variable initializer"
    initializer::Initializer = DefaultInitializer(eltype(grid))

    "Time steppers as a `NamedTuple` with `explicit` and optional `implicit` entries"
    timesteppers::Timesteppers = default_timesteppers(eltype(grid))
end

# Model interface methods

function compute_auxiliary!(state, model::SurfaceHydrologyModel)
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.canopy_interception)
    compute_auxiliary!(state, model, model.evapotranspiration)
    compute_auxiliary!(state, model, model.surface_runoff)
    return nothing
end

function compute_tendencies!(state, model::SurfaceHydrologyModel)
    compute_tendencies!(state, model, model.atmosphere)
    compute_tendencies!(state, model, model.canopy_interception)
    compute_tendencies!(state, model, model.evapotranspiration)
    compute_tendencies!(state, model, model.surface_runoff)
    return nothing
end
