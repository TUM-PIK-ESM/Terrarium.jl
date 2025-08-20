# Initializer interface

abstract type AbstractInitializer end

# Default implementations of initialize!
initialize!(state, model::AbstractModel, init::AbstractInitializer)  = nothing
initialize!(state, model::AbstractModel) = initialize!(state, model, get_initializer(model))

"""
    $SIGNATURES

Returns an `Oceananigans` `Field` initializer for the given state variable. Defaults to returning `nothing`
which will apply no initialization; note that all `Field`s are filled with zeros by default.
"""
get_field_initializer(::AbstractInitializer, ::AbstractLandGrid, ::AbstractVariable) = nothing

"""
Marker type for a no-op initializer that leaves all `Field`s set to their default values.
"""
struct DefaultInitializer <: AbstractInitializer end

"""
    $TYPEDEF

Initializer type that allows for direct specification of `Field` initializer functions.
The first argument should be the `name` of the state variable while the `init` should be a `Function`
or callable struct of the form `f(x,y,z)::Real` where `z` is the vertical coordinate, decreasing with depth.

Properties:
$TYPEDFIELDS
"""
struct VarInitializer{name, InitFunc} <: AbstractInitializer
    "Intializer function"
    init::InitFunc

    VarInitializer(name::Symbol, init) = new{name, typeof(init)}(init)
end

varname(::VarInitializer{name}) where {name} = name

"""
    $SIGNATURES

Convenience constructor of `VarInitializer` for a state variable with the given `name` and
initializer function `init`. The initializer function may optionally be declared using Julia's
`do` syntax:

```julia
var_init = VarInitializer(:var) do (x,z)
    sin(2π*z)
end
````
"""
VarInitializer(init, name::Symbol) = VarInitializer(name, init)

"""
    $SIGNATURES

Retrieves the initializer for the given variable `var` or returns `nothing` if not defined.
"""
get_field_initializer(init::VarInitializer, grid::AbstractLandGrid, var::AbstractVariable) = varname(var) == varname(init) ? init.init : nothing

"""
    $TYPEDEF

Container type that bundles multiple `AbstractInitializer`s into a single object that can be supplied to a model.
"""
struct Initializers{names, Inits<:Tuple{Vararg{AbstractInitializer}}} <: AbstractInitializer
    inits::NamedTuple{names, Inits}
end

Initializers(; inits...) = Initializers((; inits...))
Initializers(inits::VarInitializer...) = Initializers(; map(init -> varname(init) => init, inits)...)

function get_field_initializer(inits::Initializers, grid::AbstractLandGrid, var::AbstractVariable)
    field_inits = map(values(inits.inits)) do init
        get_field_initializer(init, grid, var)
    end
    # get all non-nothing values
    matched_idx = findall(!isnothing, field_inits)
    if length(matched_idx) > 1
        @warn "Found more than one matching initializer for $(varname(var)); selecting $(typeof(field_inits[matched_idx[1]]))"
    end
    return length(matched_idx) >= 1 ? field_inits[matched_idx[1]] : nothing
end
