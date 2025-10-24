struct SurfaceEnergyModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    SEB<:AbstractSurfaceEnergyBalance,
    Atmosphere<:AbstractAtmosphere,
    Initializer<:AbstractInitializer,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractSurfaceEnergyModel{NF, GridType, TimeStepper}
    "Spatial grid"
    grid::GridType

    "Atmospheric inputs"
    atmosphere::Atmosphere

    "Surface energy balance scheme"
    surface_energy_balance::SEB

    "Physical constants"
    constants::PhysicalConstants

    "State variable initializer"
    initializer::Initializer
    
    "Time stepping scheme"
    time_stepping::TimeStepper
end

function SurfaceEnergyModel(
    grid::AbstractLandGrid{NF};
    atmosphere::AbstractAtmosphere = PrescribedAtmosphere(NF),
    surface_energy_balance::AbstractSurfaceEnergyBalance = SurfaceEnergyBalance(NF),
    constants::PhysicalConstants = PhysicalConstants(NF),
    initializer::AbstractInitializer = DefaultInitializer(),
    time_stepping::AbstractTimeStepper = ForwardEuler(NF)
) where {NF}
    return SurfaceEnergyModel(grid, atmosphere, surface_energy_balance, constants, initializer, time_stepping)
end

variables(model::SurfaceEnergyModel) = tuplejoin(variables(model.atmosphere), variables(model.surface_energy_balance))

get_boundary_conditions(::SurfaceEnergyModel) = DefaultBoundaryConditions()

function compute_auxiliary!(state, model::SurfaceEnergyModel)
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.surface_energy_balance)
end

function compute_tendencies!(state, ::SurfaceEnergyModel)
    compute_tendencies!(state, model, model.atmosphere)
    compute_tendencies!(state, model, model.surface_energy_balance)
end
