# Initializer interface

abstract type AbstractInitializer end

"""
    get_field_initializer(init::AbstractInitializer, var::AbstractVariable)
"""
get_field_initializer(::AbstractInitializer, ::AbstractVariable) = nothing

"""
    $TYPEDEF

Initializer type that allows for direct specification of `Field` initializer functions.

Example:
```julia
initializer = FieldInitializer(:temperature => )
```
"""
struct FieldInitializer{VarInits<:NamedTuple} <: AbstractInitializer
    vars::VarInits
end

"""
    FieldInitializer(vars::Pair{Symbol}...)

Creates a `FieldInitializer` for the given variables. The first value of the pair is
the name of the state variable while the second value should be a `Function` or callable
struct of the form `f(x,y,z)::Real` where `z` is the vertical coordinate, decreasing with depth.

Example:
```julia
initializer = FieldInitializer(
    :temperature => (x,y,z) -> 0.01*abs(z) + exp(z)*sin(2π*z)
)
```
"""
FieldInitializer(vars::Pair{Symbol}...) = FieldInitializer((;vars...))

"""
    $SIGNATURES

Retrieves the initializer for the given variable `var` or returns `nothing` if not defined.
"""
get_initializer(init::FieldInitializer, var::AbstractVariable) = get(init.vars, varname(var), nothing)
