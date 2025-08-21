# Initializer interface

abstract type AbstractInitializer end

# Default implementations of initialize!
initialize!(state, model::AbstractModel, init::AbstractInitializer)  = nothing
initialize!(state, model::AbstractModel) = initialize!(state, model, get_initializer(model))

"""
    $SIGNATURES

Returns a named tuple of `Oceananigans` `Field` initializer functions where the keys correspond to
the names of the respective state variables.
"""
get_field_initializers(::AbstractInitializer)::NamedTuple = (;)

"""
Marker type for a no-op initializer that leaves all `Field`s set to their default values.
"""
struct DefaultInitializer <: AbstractInitializer end

"""
    $TYPEDEF

Container type that bundles multiple `AbstractInitializer`s into a single object that can be supplied to a model.
"""
struct Initializers{names, Inits, VarInits} <: AbstractInitializer
    "Component initializers"
    inits::Inits

    "Initializers for individual state variables"
    vars::NamedTuple{names, VarInits}
end

"""
Creates an `Initializers` container wrapping the given `AbstractInitializer`s `inits` as well
as any initializer functions for individual state variables in `vars`.
"""
Initializers(inits::AbstractInitializer...; vars...) = Initializers(inits, (; vars...))

function initialize!(state, model::AbstractModel, inits::Initializers{names}) where {names}
    for name in names
        set!(getproperty(state, name), inits.vars[name])
    end
    for init in inits.inits
        initialize!(state, model, init)
    end
end

get_field_initializers(inits::Initializers) = merge(init.vars, map(get_field_initializers, inits.inits)...)
