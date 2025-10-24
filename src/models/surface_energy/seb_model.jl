@kwdef struct SurfaceEnergyBalanceModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Albedo<:AbstractAlbedo,
    RadiativeFluxes<:AbstractRadiativeFluxes,
    TurbulentFluxes<:AbstractTurbulentFluxes,
    SkinTemperature<:AbstractSkinTemperature,
    Atmosphere <: AbstractAtmosphere,
    Initializer<:AbstractInitializer,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractSurfaceEnergyBalanceModel{NF, GridType, TimeStepper}
    "Spatial grid"
    grid::GridType

    "Scheme for parameterizing surface albedo"
    albedo::Albedo = ConstantAlbedo()

    "Scheme for determining the net radiation budget"
    radiative_fluxes::RadiativeFluxes = DiagnosedRadiativeFluxes()

    "Scheme for computing turbulent (sensible and latent) heat fluxes"
    turbulent_fluxes::TurbulentFluxes = DiagnosedTurbulentFluxes(eltype(grid))

    "Scheme for determining skin temperature and ground heat flux"
    skin_temperature::SkinTemperature = ImplicitSkinTemperature()

    "Atmospheric inputs"
    atmosphere::Atmosphere = PrescribedAtmosphere(grid)

    "Physical constants"
    constants::PhysicalConstants{NF} = PhysicalConstants(eltype(grid)) 

    "State variable initializer"
    initializer::Initializer = DefaultInitializer()
    
    "Time stepping scheme"
    time_stepping::TimeStepper = ImplicitSEB(eltype(grid))
end

SurfaceEnergyBalanceModel(grid::AbstractLandGrid{NF}; kwargs...) where {NF} = SurfaceEnergyBalanceModel(; grid, kwargs...)

get_boundary_conditions(::SurfaceEnergyBalanceModel) = DefaultBoundaryConditions()

variables(model::SurfaceEnergyBalanceModel) = (
    variables(model.atmosphere)...,
    variables(model.albedo)...,
    variables(model.skin_temperature)...,
    variables(model.radiative_fluxes)...,
    variables(model.turbulent_fluxes)...,
)

function compute_auxiliary!(state, model::SurfaceEnergyBalanceModel)
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.albedo)
    compute_auxiliary!(state, model, model.radiative_fluxes)
    compute_auxiliary!(state, model, model.turbulent_fluxes)
    compute_auxiliary!(state, model, model.skin_temperature)
end

function compute_tendencies!(state, model::SurfaceEnergyBalanceModel)
    compute_tendencies!(state, model, model.atmosphere)
    compute_tendencies!(state, model, model.albedo)
    compute_tendencies!(state, model, model.radiative_fluxes)
    compute_tendencies!(state, model, model.turbulent_fluxes)
    compute_tendencies!(state, model, model.skin_temperature)
end

# Specialized SEB timestepper

struct ImplicitSEB{NF,ETS<:AbstractTimeStepper{NF}} <: AbstractTimeStepper{NF}
    explicit::ETS
end

# TODO: here is a good use case for an explicit time stepper supertype
ImplicitSEB(::Type{NF}; stepper::AbstractTimeStepper=ForwardEuler(NF)) where {NF} = ImplicitSEB(stepper)

default_dt(stepper::ImplicitSEB) = default_dt(stepper.explicit)

function timestep!(
    state,
    model::SurfaceEnergyBalanceModel{NF, GR, AL, RF, TF, <:ImplicitSkinTemperature},
    stepper::ImplicitSEB,
    dt
) where {NF, GR, AL, RF, TF}
    set!(state.skin_temperature)
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    # update skin temperature according to fixed point rule
    update_skin_temperature!(state, model, model.skin_temperature)
    # update any other tendencies with nested explicit solver
    do_timestep!(state, model, stepper.explicit, dt)
end
