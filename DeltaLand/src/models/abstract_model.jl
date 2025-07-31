# AbstractModel

"""
    AbstractModel

Base type for all models.
"""
abstract type AbstractModel end

"""
    get_grid(model::AbstractModel)::AbstractGrid

Return the spatial grid associated with this `model`.
"""
function get_grid end

"""
    get_time_stepping(model::AbstractModel)::AbstractTimeStepper

Returns the time stepping scheme associated with this `model`.
"""
function get_time_stepping end

"""
    initialize(model::AbstractModel, initializer::AbstractInitializer; kwargs...)::Simulation

Initialize and return a `Simulation` based on the given `model`.
"""
function initialize end

"""
    initialize!(state, model::AbstractModel)

Fully (re-)initialize all state variables for the given `model`.
"""
function initialize! end

"""
    update_state!(state, model::AbstractModel)

Fully updates the state of `model` based on the current values of all prognostic variables in `state`.
"""
function update_state! end

"""
    compute_tendencies!(state, model::AbstractModel)

Computes tendencies for all prognostic state variables for `model` stored in the given `state`.
"""
function compute_tendencies! end

"""
    timestep!(state, model::AbstractModel, timestepper::AbstractTimeStepper, [dt = nothing])

Advance the model state by one time step, or by `dt` units of time.
"""
function timestep! end

# default method implementations

get_grid(model::AbstractModel) = model.grid

get_time_stepping(model::AbstractModel) = model.time_stepping

initialize(model::AbstractModel; kwargs...) = initialize(model, default_initializer(model); kwargs...)

initialize!(state, model::AbstractModel) = update_state!(state, model)

default_initializer(::AbstractModel) = ModelInitializer()

# TODO: define general method interfaces (as needed) for all model types
"""
    AbstractGroundModel
    
Base type for ground (e.g. soil and rock) models.
"""
abstract type AbstractGroundModel <: AbstractModel end

"""
    AbstractSoilModel

Base type for soil ground models.
"""
abstract type AbstractSoilModel <: AbstractGroundModel end

"""
    AbstractSnowModel

Base type for snow models.
"""
abstract type AbstractSnowModel <: AbstractModel end

"""
    AbstractVegetationModel

Base type for vegetation/carbon models.
"""
abstract type AbstractVegetationModel <: AbstractModel end

"""
    AbstractEnergyBalanceModel

Base type for surface energy balance models.
"""
abstract type AbstractEnergyBalanceModel <: AbstractModel end

"""
    AbstractHydrologyModel

Base type for surface hydrology models.
"""
abstract type AbstractHydrologyModel <: AbstractModel end

"""
    AbstractLandModel <: AbstractModel

Base type for full land models which couple together multiple component models.
"""
abstract type AbstractLandModel <: AbstractModel end
