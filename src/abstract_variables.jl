# Abstract variable types for declaring fields.

abstract type VarDims end

"""
    XYZ <: VarDims

Indicator type for variables that should be assigned a 3D field on their associated grid.
"""
struct XYZ <: VarDims end

"""
    XY <: VarDims

Indicator type for variables that should be assigned a 2D (lateral only) field on their associated grid.
"""
struct XY <: VarDims end

# TODO: do we need to support state variables not defined on a grid?

"""
    $TYPEDEF

Base type for state variable specification.
"""
abstract type AbstractVariable end

"""
    $TYPEDEF

Base type for prognostic variable closure relations for differential equations of the form:

```math
\\frac{\\partial g(u)}{\\partial t} = F(u)
```
where `F` represents the RHS tendency as a function of the state variable `u`, and `g(u)` is a closure or constitutive
relation that maps `u` to the physical units matching the tendency. Common examples in soil hydrothermal modeling
are temperature-enthalpy and saturation-pressure relations.
"""
abstract type AbstractClosureRelation end

"""
    varname(::AbstractVariable)
    varname(::AbstractClosureRelation)
    varname(::Pair{Symbol})

Retrieves the name of the given variable or closure. For closure relations, `varname`
should return the name of the variable returned by the closure relation.
"""
varname(var::AbstractVariable) = var.name
varname(namespace::Pair{Symbol}) = first(namespace)

"""
    vardims(::AbstractVariable)

Retrieves the grid dimensions on which this variable is defined.
"""
vardims(var::AbstractVariable) = var.dims

"""
    AuxiliaryVariable{VD<:VarDims} <: AbstractVariable

Represents an auxiliary (sometimes called "diagnostic") variable with the given `name`
and `dims` on the spatial grid.
"""
struct AuxiliaryVariable{VD<:VarDims} <: AbstractVariable
    "Name of the auxiliary variable"
    name::Symbol

    "Grid dimensions on which the variable is defined."
    dims::VD
end

struct PrognosticVariable{
    VD<:VarDims,
    TendencyVar<:AuxiliaryVariable,
    Closure<:Union{Nothing, AbstractClosureRelation}
} <: AbstractVariable
    "Name of the prognostic variable"
    name::Symbol

    "Grid dimensions on which the variable is defined."
    dims::VD

    "Closure relation for the tendency of the prognostic variable."
    closure::Closure

    "Tendency corresponding to this prognostic variable."
    tendency::TendencyVar
end

hasclosure(var::PrognosticVariable) = !isnothing(var.closure)

"""
    Namespace <: AbstractVariable

Represents a new variable namespace, typically from a subcomponent of the model.
It is (currently) assumed that tha name of the namesapce corresponds to a property
defined on the model type.
"""
struct Namespace <: AbstractVariable
    name::Symbol
end

"""
    $SIGNATURES

Convenience constructor method for `PrognosticVariable`.
"""
prognostic(name::Symbol, dims::VarDims) =
    PrognosticVariable(
        name,
        dims,
        nothing,
        AuxiliaryVariable(Symbol(name, :_tendency), dims)
    )
prognostic(name::Symbol, dims::VarDims, closure::AbstractClosureRelation) =
    PrognosticVariable(
        name,
        dims,
        # closure term
        closure,
        # tendency variable
        AuxiliaryVariable(Symbol(varname(closure), :_tendency), dims),
    )

"""
    $SIGNATURES

Convenience constructor method for `AuxiliaryVariable`.
"""
auxiliary(name::Symbol, dims::VarDims) = AuxiliaryVariable(name, dims)

"""
    $SIGNATURES

Convenience constructor method for `Namespace` provided only for consistency.
"""
namespace(name::Symbol) = Namespace(name)
