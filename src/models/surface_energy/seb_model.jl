@kwdef struct SurfaceEnergyBalanceModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Albedo<:AbstractAlbedo,
    RadiativeFluxes<:AbstractRadiativeFluxes,
    TurbulentFluxes<:AbstractTurbulentFluxes,
    SkinTemperature<:AbstractSkinTemperature,
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
    turbulent_fluxes::TurbulentFluxes = DiagnosedTurbulentFluxes()

    "Scheme for determining skin temperature and ground heat flux"
    skin_temperature::SkinTemperature = PrognosticSkinTemperature()

    "State variable initializer"
    initializer::Initializer = DefaultInitializer()
    
    "Time stepping scheme"
    time_stepping::TimeStepper = ImplicitSEB(eltype(grid))
end

SurfaceEnergyBalanceModel(grid::AbstractLandGrid{NF}; kwargs...) where {NF} = SurfaceEnergyBalanceModel(; grid, kwargs...)

function compute_auxiliary!(state, model::SurfaceEnergyBalanceModel)
    compute_auxiliary!(state, model.albedo)
    compute_auxiliary!(state, model.skin_temperature)
    compute_auxiliary!(state, model.radiative_fluxes)
    compute_auxiliary!(state, model.turbulent_fluxes)
end

function compute_tendencies!(state, model::SurfaceEnergyBalanceModel)
    compute_tendencies!(state, model.albedo)
    compute_tendencies!(state, model.skin_temperature)
    compute_tendencies!(state, model.radiative_fluxes)
    compute_tendencies!(state, model.turbulent_fluxes)
end

# Specialized SEB timestepper

struct ImplicitSEB{ETS}
    explicit::ETS
end

# TODO: here is a good use case for an explicit time stepper supertype
ImplicitSEB(::Type{NF}; stepper::AbstractTimeStepper=ForwardEuler(NF)) where {NF} = ImplicitSEB(stepper)

function timestep!(
    state,
    model::SurfaceEnergyBalanceModel{NF, GR, AL, RF, TF, <:PrognosticSkinTemperature},
    stepper::ImplicitSEB,
    dt
) where {NF, GR, AL, RF, TF}
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    # update skin temperature with custom implicit rule
    launch!(model.grid, update_skin_temperature!, state, model.grid, model.skin_temperature)
    # update any other tendencies with nested explicit solver
    do_timestep!(state, model, stepper.explicit, dt)
end
