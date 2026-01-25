"""
    $TYPEDEF

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

"""
    $TYPEDEF

Base type for sets of coupled `AbstractProcess`es. Coupled process types aggregate
together two or more sub-processes and define any 
"""
abstract type AbstractCoupledProcesses <: AbstractProcess end

"""
    $TYPEDEF

Base type for all Terrarium "models". Models are standalone representations of a system
that consist of

(i) a spatial `grid` characterizing the model domain,
(ii) zero or more `AbstractProcess`es defining the dynamics, and
(iii) an `AbstractInitializer` responsible for defining the initial state of the model.

Implementations of `AbstractModel` are required to implement, at minimum, three methods:
- [`variables`](@ref) which declares the state variables requried by the model,
- [`compute_auxiliary!`](@ref) which is responsible for computing all auxiliary (non-prognostic) variables,
- [`compute_tendencies!`](@ref) which is responsible for computing the tendencies of all prognostic variables.

Note that a default implementation of `variables` is provided which automatically collects all
variables declared by `AbstractProcess`es defined as fields (properties) of `struct`s that subtype `AbstractModel`.
"""
abstract type AbstractModel{NF, Grid<:AbstractLandGrid{NF}}  end

# Method interface for AbstractModel and AbstractProcess

"""
    variables(model::AbstractModel)

Return a `Tuple` of `AbstractVariable`s (i.e. `PrognosticVariable`, `AuxiliaryVariable`, etc.)
defined by the model or process.
"""
function variables end

"""
    initialize!(state, model::AbstractModel)

Initialize all variables defined in `state` which are defined by `model`. This defaults to simply
calling `initialize!(state, model, get_initializer(model))`.

    initialize!(state, model::AbstractModel, initializer::AbstractInitializer)

Initialize the model state variables using the corresponding `initializer`. This method only needs to be
implemented if initialization routines are necessary in addition to direct field/variable initializers.

    initialize!(state, grid, process::AbstractProcess)

Initialize all state variables associated with the given `process` defined on `model`. This method is
typically defined by the corresponding `AbstractProcess` types.
"""
function initialize! end

"""
    compute_auxiliary!(state, model::AbstractModel)

Compute updates to all auxiliary variables based on the current prognostic state of the `model`.

    compute_auxiliary!(state, grid, process:AbstractProcess, args...)

Compute all auxiliary state variables for the given `process` on `grid`.
The number and type of `args...` are permitted to vary for different process interfaces.
"""
function compute_auxiliary! end

"""
    compute_tendencies!(state, model::AbstractModel)

Compute tendencies for all prognostic state variables for `model` stored in the given `state`.
This method should be called after `compute_tendencies!`.

    compute_tendencies!(state, grid, process:AbstractProcess, args...)

Compute the tendencies of all prognostic state variables for the given `process` on `grid`.
The number and type of `args...` are permitted to vary for different process interfaces.
"""
function compute_tendencies! end

# Default method implementations

# Allow variables to be defined on any type, defaulting to an empty tuple
variables(::Any) = ()
# For AbstractCoupledProcesses and AbstractModel types, default to collecting variables on all processes contained therein
variables(obj::Union{AbstractCoupledProcesses, AbstractModel}) = mapreduce(variables, tuplejoin, processes(obj))

"""
    $TYPEDSIGNATURES

Return a tuple of `AbstractProces`es contained in the given model or coupled processes type.
Note that this is a type-stable, `@generated` function that is compiled for each argument type.
"""
@generated function processes(obj::Union{AbstractCoupledProcesses, AbstractModel})
    names = fieldnames(obj)
    types = fieldtypes(obj)
    procfields = filter(Tuple(zip(names, types))) do (name, type)
        type <: AbstractProcess
    end
    procnames = map(first, procfields)
    accessors = map(name -> :(obj.$name), procnames)
    return :(tuple($(accessors...)))
end

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
    closure!(state, model::AbstractModel)
    closure!(state, grid, process::AbstractProcess)

Apply all closure relations defined for the given `model` or `process`.
"""
closure!(state, grid, ::AbstractProcess) = nothing
function closure!(state, model::AbstractModel)
    for process in processes(model)
        closure!(state, get_grid(model), process)
    end
end

"""
    invclosure!(state, model::AbstractModel)
    invclosure!(state, grid, process::AbstractProcess)

Apply the inverse of all closure relations defined for the given `model`.
"""
invclosure!(state, grid, ::AbstractProcess) = nothing
function invclosure!(state, model::AbstractModel)
    for process in processes(model)
        invclosure!(state, get_grid(model), process)
    end
end

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
