# Initializer interface

"""
Base type for model initializers. Implementations should provide a dispatch of the `initialize!(state, model::M, init::I)` method where
`M` corresponds to the model type and `I` to the initializer. An implementation of `get_field_initializers` can also be provided which
returns a `NamedTuple` of initializer functions for individual state variable fields.
"""
abstract type AbstractInitializer end

# Default implementations of initialize!
initialize!(state, model::AbstractModel, init::AbstractInitializer) = nothing
initialize!(state, model::AbstractModel) = initialize!(state, model, get_initializer(model))

"""
"""
function initialize!(state, inits::NamedTuple{names}) where {names}
    return fastiterate(names) do name
        set!(getproperty(state, name), inits[name])
    end
end

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

Container type that bundles one or more initializer functions for individual state variable fields into
a single `AbstractInitializer`. `FieldInitializers` can also optionally wrap another `AbstractInitializer`
type whose field initializers will be merged with those in `vars`.
"""
struct FieldInitializers{names, Init <: AbstractInitializer, FieldInits} <: AbstractInitializer
    "Optional model initializer to wrap and invoke in `initialize!`"
    init::Init

    "Initializer functions for individual state variables"
    vars::NamedTuple{names, FieldInits}
end

"""
Creates a new `FieldInitializers` from the given keyword arguments where each argument corresponds
to an initializer function or value for specific state variable defined in the model. The initializers
can be any function, value, array, or `Field` that would be a valid input `x` to `Oceananigans.set!(field, x)`.
Optionally, `FieldInitializers` can wrap another 

```@example
FieldInitializers(temperature=(x,z) -> sin(2π*z), saturation_water_ice=0)
```
"""
FieldInitializers(init::AbstractInitializer = DefaultInitializer(); vars...) = FieldInitializers(init, (; vars...))

function initialize!(state, model::AbstractModel, init::FieldInitializers)
    # apply variable initializers
    initialize!(state, init.vars)
    # invoke inner initializer
    initialize!(state, model, init.init)
    return nothing
end

"""
Returns the field initializers stored in `init.vars` merged with the field initializers defined by the inner
initializer. The field initilaizers defined by the given `FieldInitializers` take precedence in the merge.
"""
get_field_initializers(init::FieldInitializers) = merge(get_field_initializers(init.init), init.vars)
