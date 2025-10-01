# AbstractModel

"""
    $TYPEDEF

Base type for all models.
"""
abstract type AbstractModel{NF, Grid<:AbstractLandGrid{NF}, TS<:AbstractTimeStepper} end

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

variables(::Any) = ()

"""
    get_grid(model::AbstractModel)::AbstractLandGrid

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

Returns the boundary conditions associated with this `model`.
"""
get_boundary_conditions(model::AbstractModel) = model.boundary_conditions

"""
    get_initializer(model::AbstractModel)::AbstractInitializer

Returns the initializer associated with this `model`.
"""
get_initializer(model::AbstractModel) = model.initializer

"""
Convenience dispatch for `timestep!` that forwards to `timestep!(state, model, get_time_stepping(model), dt)`.
"""
timestep!(state, model::AbstractModel, dt) = timestep!(state, model, get_time_stepping(model), dt)

# TODO: define general method interfaces (as needed) for all model types

"""
    $TYPEDEF
    
Base type for ground (e.g. soil and rock) models.
"""
abstract type AbstractGroundModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for soil ground models.
"""
abstract type AbstractSoilModel{NF, GR, TS} <: AbstractGroundModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for snow models.
"""
abstract type AbstractSnowModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for vegetation/carbon models.
"""
abstract type AbstractVegetationModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for surface energy balance models.
"""
abstract type AbstractEnergyBalanceModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for surface hydrology models.
"""
abstract type AbstractHydrologyModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    AbstractLandModel <: AbstractModel

Base type for full land models which couple together multiple component models.
"""
abstract type AbstractLandModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
Convenience constructor for all `AbstractLandModel` types that allows the `grid` to be passed
as the first positional argument.
"""
(::Type{Model})(grid::AbstractLandGrid; kwargs...) where {Model<:AbstractModel} = Model(; grid, kwargs...)

function Adapt.adapt_structure(to, model::AbstractModel)
    return setproperties(model, map(prop -> Adapt.adapt_structure(to, prop), getproperties(model)))
end

function Base.show(io::IO, mime::MIME"text/plain", model::AbstractModel{NF}) where {NF}
    println(io, "$(nameof(typeof(model))){$NF} with $(nameof(typeof(get_grid(model)))) on $(architecture(get_grid(model)))")
    for name in propertynames(model)
        print(io, "  $name:  ")
        # Indent property value by two spaces
        show(IOContext(io, :indent => 2), mime, getproperty(model, name))
        println(io)
    end
end
