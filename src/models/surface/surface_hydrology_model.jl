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
    SurfaceEnergyBalance<:AbstractSurfaceEnergyBalance,
    CanopyHydrology<:AbstractCanopyHydrology,
    SurfaceRunoff<:AbstractSurfaceRunoff,
    Constants<:PhysicalConstants{NF},
    Initializer<:AbstractInitializer,
} <: AbstractSurfaceHydrologyModel{NF, GridType}
    "Spatial grid type"
    grid::GridType

    "Atmospheric input configuration"
    atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Surface energy balance scheme"
    surface_energy_balance::SurfaceEnergyBalance = SurfaceEnergyBalance(eltype(grid))

    "Canopy hydrology scheme"
    canopy_hydrology::CanopyHydrology = PALADYNCanopyHydrology(eltype(grid))
    
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

get_surface_energy_balance(model::SurfaceHydrologyModel) = model.surface_energy_balance

get_canopy_hydrology(model::SurfaceHydrologyModel) = model.canopy_hydrology

get_surface_runoff(model::SurfaceHydrologyModel) = model.surface_runoff

get_constants(model::SurfaceHydrologyModel) = model.constants

# Model interface methods
variables(model::SurfaceHydrologyModel) = tuplejoin(
    variables(model.atmosphere),
    variables(model.surface_energy_balance),
    variables(model.canopy_hydrology),
    variables(model.surface_runoff),
)

processes(model::SurfaceHydrologyModel) = (
    model.atmosphere,
    model.surface_energy_balance,
    model.canopy_hydrology,
    model.surface_runoff
)

function compute_auxiliary!(state, model::SurfaceHydrologyModel)
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.surface_energy_balance)
    compute_auxiliary!(state, model, model.canopy_hydrology)
    compute_auxiliary!(state, model, model.surface_runoff)
end

function compute_tendencies!(state, model::SurfaceHydrologyModel)
    compute_tendencies!(state, model, model.atmosphere)
    compute_tendencies!(state, model, model.surface_energy_balance)
    compute_tendencies!(state, model, model.canopy_hydrology)
    compute_tendencies!(state, model, model.surface_runoff)
end

