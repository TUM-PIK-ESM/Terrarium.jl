# Initializer interface

"""
Base type for model initializers. Implementations should provide a dispatch of the `initialize!(state, model::M, init::I)` method where
`M` corresponds to the model type and `I` to the initializer. An implementation of `get_field_initializers` can also be provided which
returns a `NamedTuple` of initializer functions for individual state variable fields.
"""
abstract type AbstractInitializer{NF} end

# Default implementations of initialize!
initialize!(state, model::AbstractModel, init::AbstractInitializer) = nothing
initialize!(state, model::AbstractModel) = initialize!(state, model, get_initializer(model))
# Fallback dispatch for initialize! on process types
initialize!(state, grid, process::AbstractProcess, args...) = nothing

"""
    $TYPEDSIGNATURES

Initialize the `state` with `Field` initializers (any valid argument to `set!`) in `inits`.
"""
function initialize!(state, inits::NamedTuple{names}) where {names}
    return fastiterate(names) do name
        set!(getproperty(state, name), inits[name])
    end
end

"""
Marker type for a no-op initializer that leaves all `Field`s set to their default values.
"""
struct DefaultInitializer{NF} <: AbstractInitializer{NF} end

DefaultInitializer(::Type{NF}) where {NF} = DefaultInitializer{NF}()
