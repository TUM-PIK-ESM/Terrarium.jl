# Initializer interface

abstract type AbstractInitializer end

"""
    get_initializer(init::AbstractInitializer, var::AbstractVariable)
"""
get_initializer(::AbstractInitializer, ::AbstractVariable) = nothing

# Boundary conditions interface

abstract type AbstractBoundaryConditions end

"""
    get_boundary_conditions(bcs::AbstractBoundaryConditions, var::AbstractVariable)
"""
get_boundary_conditions(::AbstractBoundaryConditions, ::AbstractVariable) = nothing

# Field construction

function create_field(
    var::AbstractVariable,
    init::AbstractInitializer,
    bcs::AbstractBoundaryConditions,
    grid::AbstractLandGrid,
    args...;
    kwargs...
)
    FT = field_type(vardims(var))
    bcs = get_boundary_conditions(bcs, var)
    field = if !isnothing(bcs)
        FT(get_field_grid(grid), args...; boundary_conditions=bcs, kwargs...)
    else
        FT(get_field_grid(grid), args...; kwargs...)
    end
    # Apply initializer if defined
    inits = get_initializer(init, var)
    if !isnothing(inits)
        Oceananigans.set!(field, inits[name])
    end
    return field
end

field_type(::XY) = Field{Center,Center,Nothing}
field_type(::XYZ) = Field{Center,Center,Center}

"""
    FieldInitializer{VarInits<:NamedTuple}

FieldInitializer provides a simple method for bootstrapping `Field` initialization
within model implementations. This is mostly a stopgap solution and is likely to change.
"""
struct FieldInitializer{VarInits<:NamedTuple} <: AbstractInitializer
    vars::VarInits
end

"""
    FieldInitializer(vars::Pair{Symbol,<:NamedTuple}...) 

Creates a `FieldInitializer` for the given variables. The first value of the pair is
the name of the state variable while the second value should be a `NamedTuple` corresponding
of valid keyword arguments for an Oceananigans `Field`. Additionally, a property `init` may
be provided corresponding to a function which will be passed to `set!`.

Example:
```julia
initializer = FieldInitializer(
    :temperature => (x,y,z) -> sin(z)
)
```
"""
FieldInitializer(vars::Pair{Symbol}...) = FieldInitializer((;vars...))

get_initializer(init::FieldInitializer, var::AbstractVariable) = get(init.vars, varname(var), nothing)
