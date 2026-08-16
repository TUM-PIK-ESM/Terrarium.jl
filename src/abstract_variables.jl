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

function Base.summary(var::AbstractVariable)
    unitstr = varunits(var) == NoUnits ? "-" : varunits(var)
    text = "$(string(varname(var))) [$(unitstr)] on $(typeof(vardims(var)))"
    return text
end

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
Baste type for process state variables with specific intents, e.g. `prognostic`, `auxiliary`, or `input`.
"""
abstract type AbstractProcessVariable{name, VD, UT} <: AbstractVariable{name, VD, UT} end

@inline vardims(pv::AbstractProcessVariable) = vardims(pv.var)
@inline varunits(pv::AbstractProcessVariable) = varunits(pv.var)

function Base.show(io::IO, ::MIME"text/plain", var::AbstractVariable)
    units = varunits(var)
    return if units != NoUnits
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
        VD <: VarDims,
        UT <: Units,
        Var <: Variable{name, VD, UT},
        BT <: DomainSets.AbstractInterval,
        FC,
    } <: AbstractProcessVariable{name, VD, UT}
    "State variable"
    var::Var

    "Field constructor"
    ctor::FC

    "Bounds for numerical bounds of the variable"
    bounds::BT

    "Variable description"
    desc::String
end

"""
    $TYPEDEF

Represents a spatially varying input (e.g. forcing) variable with the given `name` and spatial `dims`.
Input variables can also be made to vary in time through the use of [`InputSource`](@ref)s.
"""
struct InputVariable{
        name,
        VD <: VarDims,
        UT <: Units,
        Var <: Variable{name, VD, UT},
        BT <: DomainSets.AbstractInterval,
        Def <: Union{Nothing, Number, Function},
    } <: AbstractProcessVariable{name, VD, UT}
    "State variable"
    var::Var

    "Default value or function initializer"
    default::Def

    "Variable bounds"
    bounds::BT

    "Variable description"
    desc::String
end

Adapt.adapt_structure(to, var::InputVariable) = var.var

"""
    $TYPEDEF

Represents a prognostic state variable with the given `name` and spatial `dims`. Prognostic variables
are those which are integrated by the timestepper and fully define the state of the system at any given
point in (simulation) time. From a computational perspective, they can be seen as the "roots" of the
computational graph for `update_state!`/`timestep!`. Prognostic variables generally should not be modified
by any code not belonging to the timestepper or user. They automatically define a `tendency` (auxiliary)
variable which is used to hold the value of their instantaneous time derivative computed by `compute_tendencies!`.
"""
struct PrognosticVariable{
        name,
        VD <: VarDims,
        UT <: Units,
        Var <: Variable{name, VD, UT},
        CL <: Union{Nothing, AbstractClosureRelation},
        TV <: Union{Nothing, AuxiliaryVariable},
        BT <: DomainSets.AbstractInterval,
    } <: AbstractProcessVariable{name, VD, UT}
    "State variable"
    var::Var

    "Closure relation for the tendency of the prognostic variable"
    closure::CL

    "Variable corresponding to the tendency for prognostic variables"
    tendency::TV

    "Variable bounds"
    bounds::BT

    "Variable description"
    desc::String
end

hasclosure(var::PrognosticVariable) = !isnothing(var.closure)

"""
    $TYPEDEF

Represents a new variable namespace, typically from a subcomponent of the model.
"""
struct Namespace{name, Vars}
    vars::Vars

    Namespace(name::Symbol, vars) = new{name, typeof(vars)}(vars)
end

@inline varname(ns::Namespace{name}) where {name} = varname(typeof(ns))
@inline varname(::Type{<:Namespace{name}}) where {name} = name

variables(ns::Namespace) = getfield(ns, :vars)

Base.propertynames(ns::Namespace) = (:vars, propertynames(getfield(ns, :vars))...)
Base.getproperty(ns::Namespace, name::Symbol) = name == :vars ? getfield(ns, :vars) : getproperty(getfield(ns, :vars), name)

# Variable container

"""
    $TYPEDEF

Container for abstract state variable definitions. Automatically collates and merges all variables
and namespaces passed into the constructor. Uses OrderedDicts internally to avoid NamedTuple type
explosion during initialization (each merge in a foldl creates a new distinct type), converting to
NamedTuples only at the final StateVariables construction step.
"""
struct Variables
    prognostic::OrderedDict{Symbol, AbstractVariable}
    tendencies::OrderedDict{Symbol, AuxiliaryVariable}
    auxiliary::OrderedDict{Symbol, AbstractVariable}
    inputs::OrderedDict{Symbol, AbstractVariable}
    namespaces::OrderedDict{Symbol, Namespace}
end

"""
    VarPath

Type alias for namespaced variable paths of the form `(namespace_1, ..., namespace_N, varname)`.
Used to specify the location of variables in nested namespaces.
"""
const VarPath = Tuple{Vararg{Symbol}}

"""
    varpath(name::Symbol)
    varpath(path::Pair)
    varpath(path::Tuple{Vararg{Symbol}})

Normalize the given variable name into a path of the form `(namespace_1, ..., namespace_N, varname)`.
Plain `Symbol` names correspond to variables in the root namespace, i.e. the path `(varname,)`. Namespaced
variables can be specified either as `Pair`s, e.g. `:ns1 => :ns2 => :varname`, or directly as a tuple of
`Symbol`s, e.g. `(:ns1, :ns2, :varname)`.
"""
varpath(name::Symbol) = (name,)
varpath(path::VarPath) = path
varpath(path::Pair) = (Symbol(first(path)), varpath(last(path))...)

"""
    $SIGNATURES

Wrap the given variable `var` in nested `Namespace`s according to the `path`, where `path` is the
namespace scope (i.e. the sequence of enclosing namespace names, *excluding* the variable's own name).
An empty `path` returns `var` unwrapped.
"""
with_scope(path::VarPath, var::AbstractVariable) =
    isempty(path) ? var : namespace(first(path), (with_scope(Base.tail(path), var),))

Variables(obj) = Variables(variables(obj))
Variables(vars::Variables) = vars
Variables(vars::Union{AbstractProcessVariable, Namespace}...) = Variables(vars)
function Variables(@nospecialize(vars::Tuple{Vararg{Union{AbstractProcessVariable, Namespace}}}))
    # partition variables into prognostic, auxiliary, input, and namespace groups;
    # duplicates within each group are automatically merged
    varmeta(var::AbstractVariable) = (varname(var), vardims(var), varunits(var))
    varmeta(ns::Namespace) = varname(ns)
    function register!(vardict::AbstractDict, var)
        if haskey(vardict, varname(var)) && varmeta(vardict[varname(var)]) != varmeta(var)
            error("Found incompatible duplicates of variable $(varname(var)): $(var) $(vars[varname(var)])")
        elseif !haskey(vardict, varname(var))
            vardict[varname(var)] = var
        else
            println("MERGING $(varname(var)): $var")
            vardict[varname(var)] = first(merge(vardict[varname(var)], var))
            println(vardict[varname(var)])
        end
        return nothing
    end
    # create OrderedDicts for each variable type
    prognostic_vars = OrderedDict{Symbol, AbstractVariable}()
    tendency_vars = OrderedDict{Symbol, AuxiliaryVariable}()
    auxiliary_vars = OrderedDict{Symbol, AbstractVariable}()
    input_vars = OrderedDict{Symbol, InputVariable}()
    namespaces = OrderedDict{Symbol, Namespace}()
    # register prognostic variables variables
    for var in filter(var -> isa(var, PrognosticVariable), vars)
        register!(prognostic_vars, var)
    end
    # recursively collect all closure variables; this needs to come
    # before auxiliary variable registration so that closure variable
    # Fields are available to auxiliary variable Field constructors
    for var in filter(var -> isa(var, PrognosticVariable) && hasclosure(var), vars)
        closure_vars = Variables(variables(var.closure))
        @assert isempty(closure_vars.prognostic) "Closures are not allowed to declare prognostic variables"
        @assert isempty(closure_vars.namespaces) "Closures are not allowed to declare namespaces"
        merge!(auxiliary_vars, closure_vars.auxiliary)
        merge!(input_vars, closure_vars.inputs)
    end
    # register auxiliary variables
    for var in filter(var -> isa(var, AuxiliaryVariable), vars)
        register!(auxiliary_vars, var)
    end
    # register input variables
    for var in filter(var -> isa(var, InputVariable), vars)
        name = varname(var)
        # only register input variables if they are not already declared as prognostic or auxiliary
        haskey(prognostic_vars, name) || haskey(auxiliary_vars, name) || register!(input_vars, var)
    end
    # register namespaces recursively
    for ns in filter(var -> isa(var, Namespace), vars)
        name = varname(ns)
        ns_vars = variables(ns)
        inner = Namespace(name, Variables(ns_vars))
        register!(namespaces, inner)
    end
    # register tendencies for prognostic variables
    for var in values(prognostic_vars)
        tendency_vars[varname(var)] = var.tendency
    end
    # check for illegal duplicates across variable groups
    check_duplicates(
        values(prognostic_vars)...,
        values(auxiliary_vars)...,
        values(input_vars)...,
        values(namespaces)...
    )
    return Variables(
        prognostic_vars,
        tendency_vars,
        auxiliary_vars,
        input_vars,
        namespaces,
    )
end

"""
Check for variables/namespaces with duplicate names and raise an error if duplicates are detected. Not type stable.
"""
function check_duplicates(@nospecialize(vars::Union{AbstractVariable, Namespace}...))
    names = unique(map(varname, vars))
    groups = Dict(map(n -> n => filter(==(n) ∘ varname, vars), names)...)
    for key in keys(groups)
        if length(groups[key]) > 1
            error("Found conflicting variable/namespace definitions for $key:\n$(join(groups[key], "\n"))")
        end
    end
    return
end

"""
    deduplicate_vars(vars::Tuple{Vararg{Union{AbstractVariable, Namespace}}})

Type-stable equivalent of [`deduplicate`](@ref) for tuples of `AbstractVariable`s and `Namespace`s.
"""
@generated function deduplicate_vars(vars::Tuple{Vararg{Union{AbstractVariable, Namespace}}})
    names = map(varname, vars.parameters)
    unique_idx = unique(i -> names[i], eachindex(vars.parameters))
    accessors = map(i -> :(vars[$i]), unique_idx)
    return :(tuple($(accessors...)))
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

function Base.merge(varss::Tuple{Vararg{<:Union{AbstractVariable, Namespace}}}...)
    return tuplejoin(varss...)
end

function Base.merge(vars::AbstractVariable...)
    @assert allequal(vars) "All variables must be functionally equivalent to merge"
    return first(vars)
end

function Base.merge(namespaces::Namespace...)
    names = unique(map(varname, namespaces))
    merged = map(names) do name
        group = filter(ns -> varname(ns) == name, namespaces)
        length(group) == 1 ? group[1] : Namespace(name, merge(map(ns -> Variables(variables(ns)), group)...))
    end
    return Tuple(merged)
end

function Base.propertynames(vars::Variables)
    fieldnames = (:prognostic, :tendencies, :auxiliary, :inputs, :namespaces)
    prognames = keys(getfield(vars, :prognostic))
    auxnames = keys(getfield(vars, :auxiliary))
    inputnames = keys(getfield(vars, :inputs))
    nsnames = keys(getfield(vars, :namespaces))
    return tuplejoin(fieldnames, prognames, auxnames, inputnames, nsnames)
end

function Base.getproperty(vars::Variables, name::Symbol)
    # forward getproperty calls to variable groups
    if name ∈ keys(getfield(vars, :prognostic))
        return getfield(vars, :prognostic)[name]
    elseif name ∈ keys(getfield(vars, :auxiliary))
        return getfield(vars, :auxiliary)[name]
    elseif name ∈ keys(getfield(vars, :inputs))
        return getfield(vars, :inputs)[name]
    elseif name ∈ keys(getfield(vars, :namespaces))
        return getfield(vars, :namespaces)[name]
    else
        return getfield(vars, name)
    end
end

function Base.summary(vars::Variables)
    str = "Variables(prognostic = $(keys(vars.prognostic)), auxiliary = $(keys(vars.auxiliary)), inputs = $(keys(vars.inputs)), namespaces = $(keys(vars.namespaces)))"
    return str
end

function Base.show(io::IO, vars::Variables)
    println(io, "Variables")
    println(io, "├─ Prognostic: ")
    for var in values(vars.prognostic)
        println(io, "├── $(summary(var))")
    end
    println(io, "├─ Auxiliary: ")
    for var in values(vars.auxiliary)
        println(io, "├── $(summary(var))")
    end
    println(io, "├─ Inputs: ")
    for var in values(vars.inputs)
        println(io, "├── $(summary(var))")
    end
    println(io, "├─ Namespaces:")
    for ns in values(vars.namespaces)
        println(io, "├── $(summary(ns))")
    end
    return nothing
end

# Automatically forward dispatches for `show` on tuples of variables to Variables;
# This is for the convenience of the user such that `variables(model)` pretty prints
function Base.show(
        io::IO,
        vartup::Tuple{Union{AbstractVariable, Namespace}, Vararg{Union{AbstractVariable, Namespace}}}
    )
    vars = Variables(vartup)
    show(io, vars)
    return nothing
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
@inline prognostic(name::Symbol, dims::VarDims; units = NoUnits, closure = nothing, bounds = Unbounded, desc = "") = prognostic(var(name, dims, units); closure, bounds, desc)
@inline prognostic(var::Variable; closure = nothing, bounds = Unbounded, desc = "") = PrognosticVariable(var, closure, tendency(var), bounds, desc)

"""
    $SIGNATURES

Convenience constructor method for `AuxiliaryVariable`.
"""
@inline auxiliary(name::Symbol, dims::VarDims, ctor = nothing, params = nothing; units = NoUnits, bounds = Unbounded, desc = "") = auxiliary(var(name, dims, units), ctor, params; bounds, desc)
@inline auxiliary(var::Variable, ::Nothing, ::Nothing; bounds = Unbounded, desc = "") = AuxiliaryVariable(var, nothing, bounds, desc)
@inline auxiliary(var::Variable, ctor::Function, params; bounds = Unbounded, desc = "") = AuxiliaryVariable(var, (_, grid, clock, fields) -> ctor(grid, clock, fields, params), bounds, desc)
# `KernelFunction` constructors (from `kernel`) are callable structs, not `Function`s; they define
# their own `(var, grid, clock, fields)` call convention, so store them directly as the ctor.
@inline auxiliary(var::Variable, ctor::KernelFunction, ::Nothing; bounds = Unbounded, desc = "") = AuxiliaryVariable(var, ctor, bounds, desc)

"""
    $SIGNATURES

Convenience constructor method for `InputVariable`.
"""
@inline input(name::Symbol, dims::VarDims; default = nothing, units = NoUnits, bounds = Unbounded, desc = "") = input(var(name, dims, units); default, bounds, desc)
@inline input(var::Variable; default = nothing, bounds = Unbounded, desc = "") = InputVariable(var, default, bounds, desc)

"""
    $SIGNATURES

Creates an `AuxiliaryVariable` for the tendency of a prognostic variable with the given name, dimensions, and physical units.
This constructor is primarily used internally by other constructors and does not usually need to be called by implementations of `variables`.
"""
@inline tendency(var::Variable) = auxiliary(varname(var), vardims(var), units = upreferred(varunits(var)) / u"s")

"""
    $SIGNATURES

Convenience constructor method for variable `Namespace`s.
"""
@inline namespace(name::Symbol, vars::Union{Tuple, Variables}) = Namespace(name, vars)

"""
    $SIGNATURES

Convert the given `NamedTuple` of variables into a tuple of `Namespace`s.
"""
@inline namespaces(nt::NamedTuple{names}) where {names} = map((nm, vars) -> namespace(nm, vars), names, values(nt))

"""
Alias for `Variables(vars...)`
"""
@inline variables(vars::Union{AbstractVariable, Namespace}...) = Variables(vars)

"""
Helper method that selects only prognostic variables declared on `obj`.
"""
@inline prognostic_variables(obj) = prognostic_variables(variables(obj))
@inline prognostic_variables(vars::Variables) = vars.prognostic
@inline prognostic_variables(vars::Tuple) = deduplicate_vars(Tuple(filter(var -> isa(var, PrognosticVariable), vars)))

"""
Helper method that selects only auxiliary variables declared on `obj`.
"""
@inline auxiliary_variables(obj) = auxiliary_variables(variables(obj))
@inline auxiliary_variables(vars::Variables) = vars.auxiliary
@inline auxiliary_variables(vars::Tuple) = deduplicate_vars(Tuple(filter(var -> isa(var, AuxiliaryVariable), vars)))

"""
Helper method that selects only input variables declared on `obj`.
"""
@inline input_variables(obj) = input_variables(variables(obj))
@inline input_variables(vars::Variables) = vars.inputs
@inline input_variables(vars::Tuple) = deduplicate_vars(Tuple(filter(var -> isa(var, InputVariable), vars)))

"""
Helper method that selects only closure (auxiliary) variables declared on `obj`.
"""
@inline closure_variables(obj) = closure_variables(variables(obj))
@inline function closure_variables(vars::Tuple)
    progvars = prognostic_variables(vars)
    all_closure_vars = fastmap(var -> variables(var.closure), progvars)
    return deduplicate_vars(tuplejoin(all_closure_vars...))
end

function Base.NamedTuple(vars::Tuple{Vararg{Union{AbstractVariable, Namespace}}})
    keys = map(varname, vars)
    return NamedTuple{keys}(vars)
end
