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

struct TendencyVariable{VD<:VarDims} <: AbstractVariable
    "Name of the prognostic variable"
    name::Symbol

    "Grid dimensions on which the variable is defined."
    dims::VD
end

struct PrognosticVariable{
    VD<:VarDims,
    TendencyVar<:TendencyVariable,
    Closure<:Union{Nothing,AbstractClosureRelation}
} <: AbstractVariable
    "Name of the prognostic variable"
    name::Symbol

    "Grid dimensions on which the variable is defined."
    dims::VD

    "Tendency corresponding to this prognostic variable."
    tendency::TendencyVar

    "Closure relation for the tendency of the prognostic variable."
    closure::Closure
end

hasclosure(var::PrognosticVariable) = !isnothing(var.closure)

struct AuxiliaryVariable{VD<:VarDims} <: AbstractVariable
    "Name of the auxiliary variable"
    name::Symbol

    "Grid dimensions on which the variable is defined."
    dims::VD
end

"""
    $SIGNATURES

Convenience constructor method for `PrognosticVariable`.
"""
prognostic(name::Symbol, dims::VarDims) =
    PrognosticVariable(
        name,
        dims,
        TendencyVariable(Symbol(name, :_tendency),
        dims)
    )
prognostic(name::Symbol, dims::VarDims, closure::AbstractClosureRelation) =
    PrognosticVariable(
        name,
        dims,
        TendencyVariable(Symbol(varname(closure), :_tendency)),
        closure
    )

"""
    $SIGNATURES

Convenience constructor method for `AuxiliaryVariable`.
"""
auxiliary(name::Symbol, dims::VarDims) = AuxiliaryVariable(name, dims)
