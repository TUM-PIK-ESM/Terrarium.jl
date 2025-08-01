abstract type AbstractStateVariables end

# temporary solution for holding state variables
@kwdef struct StateVariables{
    prognames, tendnames, auxnames, namespaces, cnames,
    ProgVars, TendVars, AuxVars, SubVars, Closures
} <: AbstractStateVariables
    prognostic::NamedTuple{prognames, ProgVars} = (;)
    tendencies::NamedTuple{tendnames, TendVars} = (;)
    auxiliary::NamedTuple{auxnames, AuxVars} = (;)
    namespaces::NamedTuple{namespaces, SubVars} = (;)
    closures::NamedTuple{cnames, Closures}
end

Base.propertynames(
    vars::StateVariables{prognames, tendnames, auxnames, namespaces}
) where {prognames, tendnames, auxnames, namespaces} = (
    prognames...,
    tendnames...,
    auxnames...,
    namespaces...,
    fieldnames(typeof(vars))...,
)

function Base.getproperty(
    vars::StateVariables{prognames, tendnames, auxnames},
    name::Symbol
) where {prognames, tendnames, auxnames}
    # forward getproperty calls to variable groups
    if name ∈ prognames
        return getproperty(getfield(vars, :prognostic), name)
    elseif name ∈ tendnames
        return getproperty(getfield(vars, :tendencies), name)
    elseif name ∈ auxnames
        return getproperty(getfield(vars, :auxiliary), name)
    elseif name ∈ namespaces
        return getproperty(getfield(vars, :namespaces), name)
    else
        return getfield(vars, name)
    end
end

function reset_tendencies!(state::StateVariables)
    for var in getfield(state, :tendencies)
        set!(var, zero(eltype(var)))
    end
end

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

abstract type AbstractVariable end

abstract type AbstractClosureRelation end

"""
    varname(::AbstractVariable)
    varname(::AbstractClosureRelation)

Retrieves the name of the given variable or closure. For closure relations, `varname`
should return the name of the variable returned by the closure relation.
"""
varname(var::AbstractVariable) = var.name

struct PrognosticVariable{VD<:VarDims, Closure<:Union{Nothing,AbstractClosureRelation}} <: AbstractVariable
    "Name of the prognostic variable"
    name::Symbol

    "Grid dimensions on which the variable is defined."
    dims::VD

    "Closure relation for the tendency of the prognostic variable"
    closure::Closure
end

struct AuxiliaryVariable{VD<:VarDims} <: AbstractVariable
    "Name of the auxiliary variable"
    name::Symbol

    "Grid dimensions on which the variable is defined."
    dims::VD
end

"""
Convenience constructor method for `PrognosticVariable`.
"""
prognostic(name::Symbol, dims::VarDims; closure::Union{Nothing,<:AbstractClosureRelation} = nothing) = PrognosticVariable(name, dims, closure)

"""
Convenience constructor method for `AuxiliaryVariable`.
"""
auxiliary(name::Symbol, dims::VarDims) = AuxiliaryVariable(name, dims)
