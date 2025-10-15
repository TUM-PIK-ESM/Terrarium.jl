# Abstract variable types for declaring fields.

abstract type VarDims end

"""
    XYZ <: VarDims

Indicator type for variables that should be assigned a 3D field on their associated grid.
"""
@kwdef struct XYZ{LX, LY, LZ} <: VarDims
    x::LX = Center()
    y::LY = Center()
    z::LZ = Center()
end

# Dispatch for Oceananigans `location` method
location(dims::XYZ) = (dims.x, dims.y, dims.z)

"""
    XY <: VarDims

Indicator type for variables that should be assigned a 2D (lateral only) field on their associated grid.
"""
@kwdef struct XY{LX, LY} <: VarDims
    x::LX = Center()
    y::LY = Center()
end

location(dims::XY) = (dims.x, dims.y, nothing)

# TODO: do we need to support state variables not defined on a grid?

"""
    $SIGNATURES

Infer the appropriate `VarDims` from the given `Field`.
"""
vardims(::AbstractField{LX, LY, Nothing}) where {LX, LY} = XY(LX(), LY())
vardims(::AbstractField{LX, LY, LZ}) where {LX, LY, LZ} = XYZ(LX(), LY(), LZ())

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
    getvar(::AbstractClosureRelation)

Return an `AuxiliaryVariable` corresponding to the closure variable defined by the given closure relation.
"""
function getvar end

"""
    $TYPEDEF

Base type for state variable specification.
"""
abstract type AbstractVariable{VD<:VarDims} end

"""
    $SIGNATURES

Retrieve the name of the given variable or closure. For closure relations, `varname`
should return the name of the variable returned by the closure relation.
"""
varname(var::AbstractVariable) = var.name
varname(namespace::Pair{Symbol}) = first(namespace)

"""
    $SIGNATURES

Retrieve the grid dimensions on which this variable is defined.
"""
vardims(var::AbstractVariable) = var.dims

"""
    $SIGNATURES

Retrieve the physical units for the given variable.
"""
varunits(var::AbstractVariable) = var.units

Base.:(==)(var1::AbstractVariable, var2::AbstractVariable) =
    varname(var1) == varname(var2) &&
    vardims(var1) == vardims(var2) &&
    varunits(var1) == varunits(var2)

function Base.show(io::IO, ::MIME"text/plain", var::AbstractVariable)
    units = varunits(var)
    if units != NoUnits
        println(io, "$(nameof(typeof(var))){$(typeof(vardims(var)))} with units $(string(varunits(var)))")
    else
        println(io, "$(nameof(typeof(var))){$(typeof(vardims(var)))}")
    end
end

"""
    $TYPEDEF

Represents an auxiliary (sometimes called "diagnostic") variable with the given `name`
and `dims` on the spatial grid.
"""
struct AuxiliaryVariable{VD<:VarDims, DT<:AbstractInterval, UT<:Units} <: AbstractVariable{VD}
    "Name of the auxiliary variable"
    name::Symbol

    "Grid dimensions on which the variable is defined"
    dims::VD

    "Physical untis associated with this state variable"
    units::UT

    "Interval domain on which scalar fields of this variable are defined"
    domain::DT

    "Human-readable description of this state variable"
    desc::String
end

"""
    $TYPEDEF

Represents an input (e.g. forcing) variable with the given `name` and `dims` on the spatial grid.
"""
struct InputVariable{VD<:VarDims, UT<:Units} <: AbstractVariable{VD}
    "Name of the auxiliary variable"
    name::Symbol

    "Grid dimensions on which the variable is defined"
    dims::VD

    "Physical untis associated with this state variable"
    units::UT

    "Human-readable description of this state variable"
    desc::String
end

"""
    $TYPEDEF

Represents a prognostic variable with the given `name` and `dims` on the spatial grid.
"""
struct PrognosticVariable{
    VD<:VarDims,
    UT<:Units,
    TV<:Union{Nothing, AuxiliaryVariable},
    DT<:AbstractInterval,
    CL<:Union{Nothing, AbstractClosureRelation}
} <: AbstractVariable{VD}
    "Name of the prognostic variable"
    name::Symbol

    "Grid dimensions on which the variable is defined"
    dims::VD

    "Closure relation for the tendency of the prognostic variable"
    closure::CL

    "Variable corresponding to the tendency for prognostic variables"
    tendency::TV

    "Physical untis associated with this state variable"
    units::UT

    "Interval domain on which scalar fields of this variable are defined"
    domain::DT

    "Human-readable description of this state variable"
    desc::String
end

hasclosure(var::PrognosticVariable) = !isnothing(var.closure)

"""
    $TYPEDEF

Represents a new variable namespace, typically from a subcomponent of the model.
It is (currently) assumed that tha name of the namespace corresponds to a property
defined on the model type.
"""
struct Namespace{Vars}
    name::Symbol
    vars::Vars
end

varname(ns::Namespace) = ns.name

# Variable container

"""
    $TYPEDEF

Container for abstract state variable definitions. Automatically sorts and merges all variables
and namespaces passed into the constructor.
"""
struct Variables{ProgVars, TendVars, AuxVars, InputVars, Namespaces}
    prognostic::ProgVars
    tendencies::TendVars
    auxiliary::AuxVars
    inputs::InputVars
    namespaces::Namespaces

    function Variables(vars::Union{AbstractVariable, Namespace}...)
        # partition variables into prognostic, auxiliary, input, and namespace groups;
        # duplicates within each group are automatically merged
        prognostic_vars = merge_duplicates(filter(var -> isa(var, PrognosticVariable), vars))
        auxiliary_vars = merge_duplicates(filter(var -> isa(var, AuxiliaryVariable), vars))
        input_vars = merge_duplicates(filter(var -> isa(var, InputVariable), vars))
        namespaces = merge_duplicates(filter(var -> isa(var, Namespace), vars))
        # get tendencies from prognostic variables
        tendency_vars = map(var -> var.tendency, prognostic_vars)
        # create closure variables and add to auxiliary variables
        closure_vars = map(var -> getvar(var.closure), filter(hasclosure, prognostic_vars))
        auxiliary_ext = tuplejoin(auxiliary_vars, closure_vars)
        # drop inputs with matching prognostic or auxiliary variables
        input_vars = filter(var -> var ∉ prognostic_vars && var ∉ auxiliary_vars, input_vars)
        # check for duplicates
        check_duplicates(prognostic_vars..., auxiliary_vars..., input_vars..., namespaces...)
        return new{typeof(prognostic_vars), typeof(tendency_vars), typeof(auxiliary_ext), typeof(input_vars), typeof(namespaces)}(
            prognostic_vars,
            tendency_vars,
            auxiliary_ext,
            input_vars,
            namespaces,
        )
    end
end

function check_duplicates(vars...)
    names = unique(map(varname, vars))
    groups = Dict(map(n -> n => filter(==(n) ∘ varname, vars), names)...)
    for key in keys(groups)
        if length(groups[key]) > 1
            error("Found conflicting variable/namespace definitions for $key:\n$(join(groups[key], "\n"))")
        end
    end
end

"""
    $SIGNATURES

Convenience constructors for `PrognosticVariable`.
"""
prognostic(name::Symbol, dims::VarDims, closure::Nothing=nothing; units=NoUnits, domain=RealLine(), desc="") =
    PrognosticVariable(
        name,
        dims,
        closure,
        tendency(name, dims, units),
        units,
        domain,
        desc
    )
prognostic(name::Symbol, dims::VarDims, closure::AbstractClosureRelation; units=NoUnits, domain=RealLine(), desc="") =
    PrognosticVariable(
        name,
        dims,
        # closure term
        closure,
        # tendency variable
        tendency(closure, dims),
        units,
        domain,
        desc
    )

"""
    $SIGNATURES

Convenience constructor method for `AuxiliaryVariable`.
"""
auxiliary(name::Symbol, dims::VarDims; units=NoUnits, domain=RealLine(), desc="") = AuxiliaryVariable(name, dims, units, domain, desc)

"""
    $SIGNATURES

Convenience constructor method for `InputVariable`.
"""
input(name::Symbol, dims::VarDims; units=NoUnits, desc="") = InputVariable(name, dims, units, desc)

"""
    $SIGNATURES

Creates an `AuxiliaryVariable` for the tendency of a prognostic variable with the given name, dimensions, and physical units.
This constructor is primarily used internally by other constructors and does not usually need to be called by implementations of `variables`.
"""
tendency(progname::Symbol, progdims::VarDims, progunits::Units) = AuxiliaryVariable(progname, progdims, progunits/u"s", RealLine(), "Tendency for $progname")
function tendency(closure::AbstractClosureRelation, dims::VarDims)
    # get the auxiliary closure variable
    var = getvar(closure)
    # create a tendency variable based on its name and units
    return tendency(varname(var), dims, varunits(var))
end

"""
    $SIGNATURES

Convenience constructor method for variable `Namespace`s.
"""
namespace(name::Symbol, vars::Variables) = Namespace(name, vars)
namespace(name::Symbol, vars::Tuple) = Namespace(name, Variables(vars...))

"""
Alias for `Variables(vars...)`
"""
variables(vars::Union{AbstractVariable, Namespace}...) = Variables(vars...)

"""
Type-stable and GPU-friendly placeholder for input variables in model/process `struct`s, allowing parameters/constants to
be easily replaced by input `Field`s. `Input` can be provided as the value of a struct field or in an Oceananigans BoundaryCondition type
and will result in an input variable being allocated with the given name and dimensions during model initialization.
"""
struct Input{name, units, VD, IT}
    dims::VD
    init::IT

    Input(name::Symbol, init=nothing; units::Units=NoUnits, dims::VarDims=XY()) = new{name, units, typeof(dims), typeof(init)}(dims, init)
end

# TODO: Figure out how to also include description as string; a natural workaround here would be a state variable registry
variables(in::Input{name, units}) where {name, units} = (input(name, in.dims; units),)
