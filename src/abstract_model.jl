# AbstractModel

"""
    $TYPEDEF

Base type for all models.
"""
abstract type AbstractModel{NF, Grid<:AbstractLandGrid{NF}} end

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
This method should be called after `compute_auxiliary!`.
"""
function compute_tendencies! end

# Default implementations of `AbstractModel` methods

variables(::Any) = ()

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
closure!(state, model::AbstractModel) = nothing

"""
    invclosure!(state, model::AbstractModel)

Apply the inverse of each closure relation defined for the given `model`.
"""
invclosure!(state, model::AbstractModel) = nothing

# Abstract subtypes

# TODO: define general method interfaces (as needed) for all model types

"""
    $TYPEDEF
    
Base type for ground (e.g. soil and rock) models.
"""
abstract type AbstractGroundModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for land-atmosphere energy exchange models.
"""
abstract type AbstractSurfaceEnergyModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for soil ground models.
"""
abstract type AbstractSoilModel{NF, GR} <: AbstractGroundModel{NF, GR} end

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
Convenience constructor for all `AbstractLandModel` types that allows the `grid` to be passed
as the first positional argument.
"""
(::Type{Model})(grid::AbstractLandGrid; kwargs...) where {Model<:AbstractModel} = Model(; grid, kwargs...)

function Adapt.adapt_structure(to, model::AbstractModel)
    return setproperties(model, map(prop -> Adapt.adapt_structure(to, prop), getproperties(model)))
end

function Base.show(io::IO, model::AbstractModel{NF}) where {NF}
    println(io, "$(nameof(typeof(model))){$NF} on $(architecture(get_grid(model)))")
    for name in propertynames(model)
        print(io, "├── $name:  $(summary(getproperty(model, name)))\n")
    end
end
