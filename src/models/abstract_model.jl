# AbstractModel

"""
    AbstractModel

Base type for all models.
"""
abstract type AbstractModel end

"""
    variables(model::AbstractModel)

Return a `Tuple` of `AbstractVariable`s (i.e. `PrognosticVariable`, `AuxiliaryVariable`, etc.)
defined by the model. For models that consist of one or more sub-models, variables may optionally
be grouped into namespaces by returning a `NamedTuple` where the keys correspond to the group names.
"""
function variables end

"""
    initialize!(state, model::AbstractModel)
    initialize!(state, model::AbstractModel, initializer::AbstractInitializer)

Calls `initialize!` on the `model` and its corresponding `initializer`. This method only needs to be
implemented if initialization routines are necessary in addition to direct field/variable initializers.
"""
function initialize! end

"""
    compute_auxiliary!(state, model::AbstractModel)

Computes updates to all auxiliary variables based on the current prognostic state of the `model`.
"""
function compute_auxiliary! end

"""
    compute_tendencies!(state, model::AbstractModel)

Computes tendencies for all prognostic state variables for `model` stored in the given `state`.
"""
function compute_tendencies! end

# Default getter methods for standard `AbstractModel` fields.

variables(::AbstractModel) = ()

"""
    get_grid(model::AbstractModel)::AbstractGrid

Return the spatial grid associated with this `model`.
"""
get_grid(model::AbstractModel) = model.grid

"""
    get_time_stepping(model::AbstractModel)::AbstractTimeStepper

Returns the time stepping scheme associated with this `model`.
"""
get_time_stepping(model::AbstractModel) = model.time_stepping

"""
    get_boundary_conditions(model::AbstractModel)::AbstractBoundaryConditions

Returns the time stepping scheme associated with this `model`.
"""
get_boundary_conditions(model::AbstractModel) = model.boundary_conditions

"""
    get_initializer(model::AbstractModel)::AbstractInitializer

Returns the time stepping scheme associated with this `model`.
"""
get_initializer(model::AbstractModel) = model.initializer

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

function Adapt.adapt_structure(to, model::AbstractModel)
    return setproperties(model, map(prop -> Adapt.adapt_structure(to, prop), getproperties(model)))
end
