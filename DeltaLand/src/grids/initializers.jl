abstract type AbstractInitializer end

"""
    ModelInitializer{VarInits<:NamedTuple}

ModelInitializer provides a simple method for bootstrapping `Field` initialization
within model implementations. This is mostly a stopgap solution and is likely to change.
"""
@kwdef struct ModelInitializer{VarInits<:NamedTuple} <: AbstractInitializer
    vars::VarInits = (;)
    clock::Clock = Clock(time=0.0)
end

"""
    ModelInitializer(vars::Pair{Symbol,NamedTuple}...) 

Creates a `ModelInitializer` for the given variables. The first value of the pair is
the name of the state variable while the second value should be a `NamedTuple` corresponding
of valid keyword arguments for an Oceananigans `Field`. Additionally, a property `init` may
be provided corresponding to a function which will be passed to `set!`.

Example:
```julia
initializer = ModelInitializer(
    :temperature => (boundary_conditions=ValueBoundaryCondition(0.0), init=(x,y,z) -> sin(z))
)
simulation = initialize(model, initializer)
```
"""
ModelInitializer(vars::Pair{Symbol,NamedTuple}...; clock=Clock(time=0.0)) = ModelInitializer(; vars=(;vars...), clock)

function initialize(
    ::Type{T},
    init::ModelInitializer,
    grid::AbstractLandGrid,
    name::Symbol,
    args...;
    kwargs...
)
    # get initialization fields if present
    fieldinfo = get(init.inits, name, (;))
    # filter out reserved 'init' field and interpret the rest as Field kwargs
    field_kwargs = Base.structdiff(fieldinfo, NamedTuple{(:init,)})
    field = T(grid, args...; field_kwargs..., kwargs...)
    :init ∈ fieldinfo && Oceananigans.set!(field, fieldinfo.init)
    return field
end

# TODO: also need to handle non-centered fields (e.g. boundary conditions)
initialize1D(init::ModelInitializer, grid::AbstractLandGrid, name::Symbol, args...; kwargs...) = initialize(Field{Nothing,Nothing,Center}, init, grid, args...; kwargs...)
initialize2D(init::ModelInitializer, grid::AbstractLandGrid, name::Symbol, args...; kwargs...) = initialize(Field{Center,Center,Nothing}, init, grid, args...; kwargs...)
initialize3D(init::ModelInitializer, grid::AbstractLandGrid, name::Symbol, args...; kwargs...) = initialize(CenterField, grid, name, args...; kwargs...)

