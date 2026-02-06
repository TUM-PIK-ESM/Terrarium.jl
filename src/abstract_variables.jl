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
Oceananigans.location(dims::XYZ) = (dims.x, dims.y, dims.z)

"""
    XY <: VarDims

Indicator type for variables that should be assigned a 2D (lateral only) field on their associated grid.
"""
@kwdef struct XY{LX, LY} <: VarDims
    x::LX = Center()
    y::LY = Center()
end

Oceananigans.location(dims::XY) = (dims.x, dims.y, nothing)

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
Base type for state variable placeholder types.
"""
abstract type AbstractVariable{name, VD, UT} end

"""
    $SIGNATURES

Retrieve the name of the given variable or closure. For closure relations, `varname`
should return the name of the variable returned by the closure relation.
"""
@inline varname(::AbstractVariable{name}) where {name} = name
@inline varname(::Type{<:AbstractVariable{name}}) where {name} = name
@inline varname(namespace::Pair{Symbol}) = first(namespace)

"""
    $SIGNATURES

Retrieve the grid dimensions on which this variable is defined.
"""
@inline vardims(var::AbstractVariable) = var.dims
@inline vardims(::Type{<:AbstractVariable{name, VD}}) where {name, VD} = VD

"""
    $SIGNATURES

Retrieve the physical units for the given variable.
"""
@inline varunits(var::AbstractVariable) = var.units
@inline varunits(::Type{<:AbstractVariable{name, VD, UT}}) where {name, VD, UT} = UT

# Test equality between variables by their names, dimensions, and physical units
Base.:(==)(var1::AbstractVariable, var2::AbstractVariable) =
    varname(var1) == varname(var2) &&
    vardims(var1) == vardims(var2) &&
    varunits(var1) == varunits(var2)

"""
    $TYPEDEF

Represents metadata for a generic state variable with the given `name` and spatial `dims`.
"""
struct Variable{name, VD, UT} <: AbstractVariable{name, VD, UT}
    "Variable dimensions"
    dims::VD

    "Physical units"
    units::UT

    Variable(name::Symbol, dims::VarDims, units::Units = NoUnits) = new{name, typeof(dims), typeof(units)}(dims, units)
end

"""
Baste type for process state variables with specific intents, e.g. `prognostic`, `auxiliary`, or `input`.
"""
abstract type AbstractProcessVariable{name, VD, UT} <: AbstractVariable{name, VD, UT} end

@inline vardims(pv::AbstractProcessVariable) = vardims(pv.var)
@inline varunits(pv::AbstractProcessVariable) = varunits(pv.var)

function Base.show(io::IO, ::MIME"text/plain", var::AbstractVariable)
    units = varunits(var)
    if units != NoUnits
        println(io, "$(nameof(typeof(var))) $(varname(var)) with dimensions $(typeof(vardims(var))) and units $(string(varunits(var)))")
    else
        println(io, "$(nameof(typeof(var))) $(varname(var)) with dimensions $(typeof(vardims(var)))")
    end
end

"""
    $TYPEDEF

Represents an auxiliary (a.k.a "diagnostic") state variable with the given `name`
and spatial `dims`. Auxiliary variables are those which are diagnosed directly or
indirectly from the values of one or more prognostic variables.
"""
struct AuxiliaryVariable{
    name,
    VD<:VarDims,
    UT<:Units,
    Var<:Variable{name, VD, UT},
    DT<:AbstractInterval,
    FC<:Union{Nothing, Function}
} <: AbstractProcessVariable{name, VD, UT}
    "State variable"
    var::Var

    "Field constructor"
    ctor::FC

    "Variable domain"
    domain::DT

    "Variable description"
    desc::String
end

"""
    $TYPEDEF

Represents an input (e.g. forcing) variable with the given `name` and spatial `dims`.
"""
struct InputVariable{
    name,
    VD<:VarDims,
    UT<:Units,
    Var<:Variable{name, VD, UT},
    DT<:AbstractInterval,
    Def<:Union{Nothing, Number, Function}
} <: AbstractProcessVariable{name, VD, UT}
    "State variable"
    var::Var

    "Default value or function initializer"
    default::Def

    "Variable domain"
    domain::DT

    "Variable description"
    desc::String
end

"""
    $TYPEDEF

Represents a prognostic state variable with the given `name` and spatial `dims`.
"""
struct PrognosticVariable{
    name,
    VD<:VarDims,
    UT<:Units,
    Var<:Variable{name, VD, UT},
    CL<:Union{Nothing, AbstractClosureRelation},
    TV<:Union{Nothing, AuxiliaryVariable},
    DT<:AbstractInterval
} <: AbstractProcessVariable{name, VD, UT}
    "State variable"
    var::Var

    "Closure relation for the tendency of the prognostic variable"
    closure::CL

    "Variable corresponding to the tendency for prognostic variables"
    tendency::TV

    "Variable domain"
    domain::DT

    "Variable description"
    desc::String
end

hasclosure(var::PrognosticVariable) = !isnothing(var.closure)

# Variable container

"""
    $TYPEDEF

Container for abstract state variable definitions. Automatically collates and merges all variables
and namespaces passed into the constructor.
"""
struct Variables{ProgVars, TendVars, AuxVars, InputVars, Namespaces}
    prognostic::ProgVars
    tendencies::TendVars
    auxiliary::AuxVars
    inputs::InputVars
    namespaces::Namespaces
end

"""
    $TYPEDEF

Represents a new variable namespace, typically from a subcomponent of the model.
"""
struct Namespace{name, Vars}
    vars::Vars

    Namespace(name::Symbol, vars::Variables) = new{name, typeof(vars)}(vars)
end

@inline varname(::Namespace{name}) where {name} = name

Variables(obj) = Variables(variables(obj))
Variables(vars::Union{AbstractProcessVariable, Namespace}...) = Variables(vars)
function Variables(@nospecialize(vars::Tuple{Vararg{Union{AbstractProcessVariable, Namespace}}}))
    # partition variables into prognostic, auxiliary, input, and namespace groups;
    # duplicates within each group are automatically merged
    varinfo(var::AbstractVariable) = (varname(var), vardims(var), varunits(var))
    varinfo(ns::Namespace) = varname(ns)
    # Note that we intentionally use `deduplicate` here instead of `deduplicate_vars` to
    # avoid unnecessary compilation overhead; this is performance non-critical code anyway.
    prognostic_vars = deduplicate(varinfo, filter(var -> isa(var, PrognosticVariable), vars))
    auxiliary_vars = deduplicate(varinfo, filter(var -> isa(var, AuxiliaryVariable), vars))
    input_vars = deduplicate(varinfo, filter(var -> isa(var, InputVariable), vars))
    namespaces = deduplicate(varinfo, filter(var -> isa(var, Namespace), vars))
    # get tendencies from prognostic variables
    tendency_vars = map(var -> var.tendency, prognostic_vars)
    # create closure variables and prepend to tuple of auxiliary variables;
    # note that the order matters here since Field constructors will be called in the order
    # that they appear in the var tuples.
    closure_vars = mapreduce(var -> variables(var.closure), tuplejoin, filter(hasclosure, prognostic_vars), init=())
    auxiliary_vars = deduplicate(varinfo, tuplejoin(closure_vars, auxiliary_vars))
    # drop inputs with matching prognostic or auxiliary variables
    input_vars = filter(var -> var ∉ prognostic_vars && var ∉ auxiliary_vars, input_vars)
    # check for duplicates
    check_duplicates(prognostic_vars..., auxiliary_vars..., input_vars..., namespaces...)
    # create named tuples for each variable type
    prognostic_nt = (; map(var -> varname(var) => var, prognostic_vars)...)
    tendency_nt = (; map(var -> varname(var) => var, tendency_vars)...)
    auxiliary_nt = (; map(var -> varname(var) => var, auxiliary_vars)...)
    input_nt = (; map(var -> varname(var) => var, input_vars)...)
    namespace_nt = (; map(var -> varname(var) => var, namespaces)...)
    return Variables(
        prognostic_nt,
        tendency_nt,
        auxiliary_nt,
        input_nt,
        namespace_nt,
    )
end

"""
    $TYPEDSIGNATURES

Removes duplicate variables from the given tuple of `vars` by name.
"""
@generated function deduplicate_vars(vars::Tuple{Vararg{AbstractVariable}})
    names = map(varname, vars.parameters)
    unique_idx = unique(i -> names[i], eachindex(vars.parameters))
    accessors = map(i -> :(vars[$i]), unique_idx)
    return :(tuple($(accessors...)))
end

"""
Check for variables/namespaces with duplicate names and raise an error if duplicates are detected.
"""
function check_duplicates(vars::Union{AbstractVariable, Namespace}...)
    names = unique(map(varname, vars))
    groups = Dict(map(n -> n => filter(==(n) ∘ varname, vars), names)...)
    for key in keys(groups)
        if length(groups[key]) > 1
            error("Found conflicting variable/namespace definitions for $key:\n$(join(groups[key], "\n"))")
        end
    end
end

"""
Merges all of the given `Variables` containers into a single container.
"""
function Base.merge(varss::Variables...)
    allvars = map(varss) do vars
        tuplejoin(
            values(vars.prognostic),
            values(vars.auxiliary),
            values(vars.inputs),
            values(vars.namespaces)
        )
    end
    return Variables(reduce(tuplejoin, allvars))
end

"""
    $SIGNATURES

Convenience constructor for `Variable`.
"""
@inline var(name::Symbol, dims::VarDims, units::Units = NoUnits) = Variable(name, dims, units)

"""
    $SIGNATURES

Convenience constructors for `PrognosticVariable`.
"""
@inline prognostic(name::Symbol, dims::VarDims; units=NoUnits, closure=nothing, domain=RealLine(), desc="") = prognostic(var(name, dims, units); closure, domain, desc)
@inline prognostic(var::Variable; closure=nothing, domain=RealLine(), desc="") = PrognosticVariable(var, closure, tendency(var), domain, desc)

"""
    $SIGNATURES

Convenience constructor method for `AuxiliaryVariable`.
"""
@inline auxiliary(name::Symbol, dims::VarDims, ctor = nothing, params = nothing; units = NoUnits, domain = RealLine(), desc = "") = auxiliary(var(name, dims, units), ctor, params; domain, desc)
@inline auxiliary(var::Variable, ::Nothing, ::Nothing; domain=RealLine(), desc="") = AuxiliaryVariable(var, nothing, domain, desc)
@inline auxiliary(var::Variable, ctor::Function, params; domain=RealLine(), desc="") = AuxiliaryVariable(var, (args...) -> ctor(params, args...), domain, desc)

"""
    $SIGNATURES

Convenience constructor method for `InputVariable`.
"""
@inline input(name::Symbol, dims::VarDims; default = nothing, units = NoUnits, domain = RealLine(), desc="") = input(var(name, dims, units); default, domain, desc)
@inline input(var::Variable; default = nothing, domain = RealLine(), desc="") = InputVariable(var, default, domain, desc)

"""
    $SIGNATURES

Creates an `AuxiliaryVariable` for the tendency of a prognostic variable with the given name, dimensions, and physical units.
This constructor is primarily used internally by other constructors and does not usually need to be called by implementations of `variables`.
"""
@inline tendency(var::Variable) = auxiliary(varname(var), vardims(var), units=upreferred(varunits(var))/u"s")

"""
    $SIGNATURES

Convenience constructor method for variable `Namespace`s.
"""
@inline namespace(name::Symbol, vars::Variables) = Namespace(name, vars)
@inline namespace(name::Symbol, vars::Tuple) = Namespace(name, Variables(vars...))

"""
Alias for `Variables(vars...)`
"""
@inline variables(vars::Union{AbstractVariable, Namespace}...) = Variables(vars)

"""
Helper method that selects only prognostic variables declared on `obj`.
"""
@inline prognostic_variables(obj) = prognostic_variables(variables(obj))
@inline prognostic_variables(vars::Tuple) = deduplicate_vars(Tuple(filter(var -> isa(var, PrognosticVariable), vars)))

"""
Helper method that selects only auxiliary variables declared on `obj`.
"""
@inline auxiliary_variables(obj) = auxiliary_variables(variables(obj))
@inline auxiliary_variables(vars::Tuple) = deduplicate_vars(Tuple(filter(var -> isa(var, AuxiliaryVariable), vars)))

"""
Helper method that selects only input variables declared on `obj`.
"""
@inline input_variables(obj) = input_variables(variables(obj))
@inline input_variables(vars::Tuple) = deduplicate_vars(Tuple(filter(var -> isa(var, InputVariable), vars)))

"""
Helper method that selects only closure (auxiliary) variables declared on `obj`.
"""
@inline closure_variables(obj) = closure_variables(variables(obj))
@inline function closure_variables(vars::Tuple)
    progvars = prognostic_variables(vars)
    closure_vars = mapreduce(var -> variables(var.closure), tuplejoin, progvars, init=())
    return deduplicate_vars(closure_vars)
end


function Base.NamedTuple(vars::Tuple{Vararg{Union{AbstractVariable, Namespace}}})
    keys = map(varname, vars)
    return NamedTuple{keys}(vars)
end
