"""
    AbstractProcess

Base type for all "processes". Implementations of `AbstractProcess` define equations,
state variables, and parameterizations which characterize the dynamics of a system at
for any given transient state. Note that processes should be largely agnostic to the
details regarding spatial and temporal discretization of the model; i.e. they should
not require specification of a specific grid or time stepping scheme but rather should
be able to operate on any given set of `Field`s and parameters representing the state of
a model at any point in time. Note that process types may also wrap/orchestrate one or
more other process types.
"""
abstract type AbstractProcess end

# AbstractModel interface

"""
    $TYPEDEF

Base type for all Terrarium "models". Models are collections of one or more processes
that additionally define

(i) a spatial `grid` characterizing the model domain, and
(ii) an `AbstractInitializer` responsible for defining the initial state of the model.
"""
abstract type AbstractModel{NF, Grid<:AbstractLandGrid{NF}} <: AbstractProcess end

"""
    variables(model::AbstractModel)
    variables(model::AbstractModel)

Return a `Tuple` of `AbstractVariable`s (i.e. `PrognosticVariable`, `AuxiliaryVariable`, etc.)
defined by the model or process.
"""
function variables end

"""
    get_processes(model::AbstractModel)

Return a tuple of `AbstractProcess` types defiend by this `model`.
"""
function get_processes end

"""
    initialize!(state, model::AbstractModel)

Initialize all variables defined in `state` which are defined by `model`. This defaults to simply
calling `initialize!(state, model, get_initializer(model))`.

    initialize!(state, model::AbstractModel, initializer::AbstractInitializer)

Initialize the model state variables using the corresponding `initializer`. This method only needs to be
implemented if initialization routines are necessary in addition to direct field/variable initializers.

    initialize!(state, model::AbstractModel, process::AbstractProcess)

Initialize all state variables associated with the given `process` defined on `model`. This method is
typically defined by the corresponding `AbstractProcess` types.
"""
function initialize! end

"""
    compute_auxiliary!(state, model::AbstractModel)

Compute updates to all auxiliary variables based on the current prognostic state of the `model`.

    compute_auxiliary!(state, model::AbstractModel, process:AbstractProcess)

Compute updates to auxiliary variables for the given `process` defined on `model`.
"""
function compute_auxiliary! end

"""
    compute_tendencies!(state, model::AbstractModel)

Compute tendencies for all prognostic state variables for `model` stored in the given `state`.
This method should be called after `compute_tendencies!`.

    compute_tendencies!(state, model::AbstractModel, process:AbstractProcess)

Compute tendencies of all prognostic state variables for the given `process` defined on `model`.
"""
function compute_tendencies! end

# Default implementations of `AbstractModel` methods

# allow variables to be defined on any type
variables(::Any) = ()

"""
    get_processes(::AbstractModel)

Return a tuple of `AbstractProcess`es defind by the given `model`.
"""
get_processes(::AbstractModel) = ()

"""
    get_grid(model::AbstractModel)::AbstractLandGrid

Return the spatial grid associated with this `model`.
"""
get_grid(model::AbstractModel) = model.grid

"""
    get_initializer(model::AbstractModel)::AbstractInitializer

Returns the initializer associated with this `model`.
"""
get_initializer(model::AbstractModel) = model.initializer

"""
    get_closures(model::AbstractModel)

Return all closure relations defined for the given `model`.
"""
get_closures(model::AbstractModel) = ()

"""
    closure!(state, model::AbstractModel)

Apply each closure relation defined for the given `model`.
"""
function closure!(state, model::AbstractModel)
    for closure in get_closures(model)
        closure!(state, model, closure)
    end
end

"""
    invclosure!(state, model::AbstractModel)

Apply the inverse of each closure relation defined for the given `model`.
"""
function invclosure!(state, model::AbstractModel)
    for closure in get_closures(model)
        invclosure!(state, model, closure)
    end
end

# Default implementation for processes, also allowing for dispatches on `nothing`
# TODO: Is this a good idea? Should we force users to *always* define these methods?
initialize!(state, model, ::Union{Nothing, AbstractProcess}) = nothing
compute_auxiliary!(state, model, ::Union{Nothing, AbstractProcess}) = nothing
compute_tendencies!(state, model, ::Union{Nothing, AbstractProcess}) = nothing
closure!(state, model, ::Union{Nothing, AbstractProcess}) = nothing
invclosure!(state, model, ::Union{Nothing, AbstractProcess}) = nothing

# AbstractModel subtypes

# TODO: define general method interfaces (as needed) for all model types

"""
    $TYPEDEF
    
Base type for ground (e.g. soil and rock) models.
"""
abstract type AbstractGroundModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for soil ground models.
"""
abstract type AbstractSoilModel{NF, GR} <: AbstractGroundModel{NF, GR} end

"""
    $TYPEDEF

Base type for land-atmosphere energy exchange models.
"""
abstract type AbstractSurfaceEnergyModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for snow models.
"""
abstract type AbstractSnowModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for vegetation/carbon models.
"""
abstract type AbstractVegetationModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for surface hydrology models.
"""
abstract type AbstractHydrologyModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    AbstractLandModel <: AbstractModel

Base type for full land models which couple together multiple component models.
"""
abstract type AbstractLandModel{NF, GR} <: AbstractModel{NF, GR} end

"""
Convenience constructor for all `AbstractModel` types that allows the `grid` to be passed
as the first positional argument.
"""
(::Type{Model})(grid::AbstractLandGrid, args...; kwargs...) where {Model<:AbstractModel} = Model(args...; grid, kwargs...)

function Base.show(io::IO, model::AbstractModel{NF}) where {NF}
    println(io, "$(nameof(typeof(model))){$NF} on $(architecture(get_grid(model)))")
    for name in propertynames(model)
        print(io, "├── $name:  $(summary(getproperty(model, name)))\n")
    end
end
