abstract type AbstractStateVariables end

"""
    $TYPEDEF

Container type for all `Field`s corresponding to state variables defined by a model.
`StateVariables` partitions the fields into three categories: prognostic, tendencies, and
auxiliary. Prognostic variables are those which characterize the state of the system and
are assigned tendencies to be integrated by the timestepper. Auxiliary fields are additional
state variables derived from the prognostic state variables but which are conditionally
independent of their values at the previous time step given the current prognostic state.
It is worth noting that tendencies are also treated internally as auxiliary variables;
however, they are assigned their own category here since they need to be handled separately
by the timestepping scheme.
"""
struct StateVariables{
    NF,
    prognames, tendnames, auxnames, inputnames, nsnames,
    ProgFields, TendFields, AuxFields, InputFields, Namespaces,
    ClockType
} <: AbstractStateVariables
    prognostic::NamedTuple{prognames, ProgFields}
    tendencies::NamedTuple{tendnames, TendFields}
    auxiliary::NamedTuple{auxnames, AuxFields}
    inputs::NamedTuple{inputnames, InputFields}
    namespaces::NamedTuple{nsnames, Namespaces}
    clock::ClockType

    function StateVariables(
        ::Type{NF},
        prognostic::NamedTuple{prognames, ProgFields},
        tendencies::NamedTuple{tendnames, TendFields},
        auxiliary::NamedTuple{auxnames, AuxFields},
        inputs::NamedTuple{inputnames, InputFields},
        namespaces::NamedTuple{nsnames, Namespaces},
        clock::ClockType
    ) where {NF, prognames, tendnames, auxnames, inputnames, nsnames,
             ProgFields, TendFields, AuxFields, InputFields, Namespaces, ClockType}
        return new{NF, prognames, tendnames, auxnames, inputnames, nsnames,
                   ProgFields, TendFields, AuxFields, InputFields, Namespaces, ClockType}(
            prognostic,
            tendencies,
            auxiliary,
            inputs,
            namespaces,
            clock
        )
    end
end

# Allow reconstruction from properties
ConstructionBase.constructorof(::Type{StateVariables{NF}}) where {NF} = (args...) -> StateVariables(NF, args...)

"""
    update_state!(state::StateVariables, model::AbstractModel, inputs::InputSources; compute_tendencies = true)

Update the `state` for the given `model` and `inputs`; this includes calling `update_inputs!` and
`fill_halo_regions!` followed by `compute_auxiliary!` and `compute_tendencies!`, if `compute_tendencies = true`.
"""
function Oceananigans.TimeSteppers.update_state!(state::StateVariables, model::AbstractModel, inputs::InputSources; compute_tendencies = true)
    reset_tendencies!(state)
    update_inputs!(state, inputs)
    fill_halo_regions!(state)
    compute_auxiliary!(state, model)
    if compute_tendencies
        compute_tendencies!(state, model)
    end
end

"""
Invoke `fill_halo_regions!` for all fields in `state`.
"""
function Oceananigans.BoundaryConditions.fill_halo_regions!(state::StateVariables)
    # fill_halo_regions! for all prognostic variables
    for field in state.prognostic
        fill_halo_regions!(field, state.clock, state)
    end
    # fill_halo_regions! for all auxiliary variables
    for field in state.auxiliary
        fill_halo_regions!(field, state.clock, state)
    end
    # recurse over namespaces
    for ns in state.namespaces
        fill_halo_regions!(ns)
    end
end

"""
Update input variables from the given input `sources`.
"""
function update_inputs!(state::StateVariables, sources::InputSources)
    # update inputs in current namespace
    update_inputs!(state.inputs, sources, state.clock)
    # recursively update namespaces
    for ns in state.namespaces
        update_inputs!(ns, sources)
    end
end

"""
Reset all tendencies in `state` to zero.
"""
function reset_tendencies!(state::StateVariables)
    # reset all tendency fields
    fastiterate(state.tendencies) do field
        set!(field, zero(eltype(field)))
    end
    # recurse over namespaces
    fastiterate(state.namespaces) do ns
        reset_tendencies!(ns)
    end
end

"""
    get_fields(state, queries::Union{Symbol, Pair}...)

Retrieves fields with names given in `queries` and returns them in a `NamedTuple`. Each argument
in `queries` can either be a `Symbol` corresponding to a field/variable defined in the namespace
of `state` or a `Pair{Symbol, Tuple}` where the key is the child namespace and the value is a
tuple of queries from that namespace.

!!! warning
    This method relies on runtime dispatch and thus should not be used in performance-critical code.
    If you need to query fields for specific sets of variables or components, use one of the
    type-stable variants instead.

```julia
# initialize model
state = initialize(model)
# get the temperature and saturation_water_ice fields
fields = get_fields(state, :temperature, :saturation_water_ice)
# extract temperature as well as variables from a namespace
nested_fields = get_fields(state, :temperature, :namespace => (:subvar1, :subvar2))
```
"""
function get_fields(state, queries::Union{Symbol, Pair}...)
    fields = map(queries) do query
        if isa(query, Symbol)
            query => getproperty(state, query)
        else isa(query, Pair)
            key, value = query
            @assert isa(value, Tuple) "namespaces queries must be given as tuples"
            key => get_fields(getproperty(state, key), value...)
        end
    end
    return (; fields...)
end

"""
    $TYPEDSIGNATURES

Retrieves the `Field` from `state` matching the `name` of the given variable.
"""
@inline get_field(state, ::AbstractVariable{name}) where {name} = getproperty(state, name)

"""
    $TYPEDSIGNATURES

Retrieves all `Field`s from `state` matching the names of the given variables.
"""
@inline function get_fields(state, vars::Tuple{Vararg{AbstractVariable}})
    names = map(varname, vars)
    matched_fields = fastmap(vars) do var
        get_field(state, var)
    end
    return NamedTuple{names}(matched_fields)
end

"""
    $TYPEDSIGNATURES

Retrieves all `Field`s from `state` corresponding to prognostic variables defined on the given `components`.
"""
@inline function prognostic_fields(state, components::Union{AbstractModel, AbstractProcess}...)
    return get_fields(state, mapreduce(prognostic_variables, tuplejoin, components))
end

"""
    $TYPEDSIGNATURES

Retrieves all `Field`s from `state` corresponding to auxiliary variables defined on the given `components`.
"""
@inline function auxiliary_fields(state, components::Union{AbstractModel, AbstractProcess}...)
    return get_fields(state, mapreduce(auxiliary_variables, tuplejoin, components))
end

"""
    $TYPEDSIGNATURES

Retrieves all `Field`s from `state` corresponding to input variables defined on the given `components`.
"""
@inline function input_fields(state, components::Union{AbstractModel, AbstractProcess}...)
    return get_fields(state, mapreduce(input_variables, tuplejoin, components))
end

# Initialization of StateVariables from models and processes

"""
    initialize(
        process::AbstractProcess,
        grid::AbstractLandGrid{NF};
        clock = Clock(time=zero(NF)),
        input_variables = (),
        boundary_conditions = (;),
        fields = (;)
    ) where {NF}

Initialize a `StateVariables` data structure containing `Field`s for all variables defined by `model`,
initialized on its associated `grid`. Any predefined `boundary_conditions` and `fields` will be passed
through to `initialize` for each variable.
"""
function initialize(
    model::AbstractModel{NF};
    clock = Clock(time=zero(NF)),
    input_variables = (),
    boundary_conditions = (;),
    fields = (;)
) where {NF}
    vars = Variables(tuplejoin(variables(model), input_variables))
    state = initialize(vars, model.grid; clock, boundary_conditions, fields)
    return state
end

"""
    initialize(
        process::AbstractProcess,
        grid::AbstractLandGrid{NF};
        clock = Clock(time=zero(NF)),
        input_variables = (),
        boundary_conditions = (;),
        fields = (;)
    ) where {NF}

Initialize a `StateVariables` data structure containing `Field`s defined on the given `grid`
for all variables defined by `process`. Any predefined `boundary_conditions` and `fields` will be passed through to `initialize`
for each variable.
"""
function initialize(
    process::AbstractProcess,
    grid::AbstractLandGrid{NF};
    clock = Clock(time=zero(NF)),
    input_variables = (),
    boundary_conditions = (;),
    fields = (;)
) where {NF}
    vars = Variables(tuplejoin(variables(process), input_variables))
    state = initialize(vars, grid; clock, boundary_conditions, fields)
    return state
end

# Initialization from variable metadata

"""
    initialize(
        vars::Variables,
        grid::AbstractLandGrid{NF};
        clock::Clock = Clock(time=0.0),
        boundary_conditions = (;),
        fields = (;)
    ) where {NF}

Initialize a `StateVariables` data structure containing `Field`s defined on the given `grid`
for all variables in `vars`. Any predefined `boundary_conditions` and `fields` will be passed through to `initialize`
for each variable.
"""
function initialize(
    vars::Variables,
    grid::AbstractLandGrid{NF};
    clock::Clock = Clock(time=0.0),
    boundary_conditions = (;),
    fields = (;)
) where {NF}
    # Initialize Fields for each variable group, if they are not already given in the user defined `fields`
    input_fields = initialize(vars.inputs, grid, clock, boundary_conditions, fields)
    tendency_fields = initialize(vars.tendencies, grid, clock, boundary_conditions, fields)
    prognostic_fields = initialize(vars.prognostic, grid, clock, boundary_conditions, merge(fields, input_fields))
    auxiliary_fields = initialize(vars.auxiliary, grid, clock, boundary_conditions, merge(fields, input_fields, prognostic_fields))
    # recursively initialize state variables for each namespace
    namespaces = map(vars.namespaces) do ns
        ns_bcs = get(boundary_conditions, varname(ns), (;))
        ns_fields = get(fields, varname(ns), (;))
        initialize(ns.vars, grid; clock, boundary_conditions=ns_bcs, fields=ns_fields)
    end
    # construct and return StateVariables
    return StateVariables(
        NF,
        prognostic_fields,
        tendency_fields,
        auxiliary_fields,
        input_fields,
        namespaces,
        clock,
    )
end

# Base case: empty named tuples
initialize(::NamedTuple{(), Tuple{}}, ::AbstractLandGrid, ::Clock, ::NamedTuple, ::NamedTuple) = (;)

"""
    initialize(
        vars::NamedTuple{names, <:Tuple{Vararg{AbstractVariable}}},
        grid::AbstractLandGrid,
        clock::Clock,
        boundary_conditions::NamedTuple,
        fields::NamedTuple
    ) where {names}

Initialize `Field`s on `grid` for each of the variables in the given named tuple `vars`.
Any predefined `boundary_conditions` and `fields` will be passed through to `initialize`
for each variable.
"""
function initialize(
    vars::NamedTuple{names, <:Tuple{Vararg{AbstractVariable}}},
    grid::AbstractLandGrid,
    clock::Clock,
    boundary_conditions::NamedTuple,
    fields::NamedTuple
) where {names}
    # Initialize or retrieve Fields for each variable in `var`, accumulating the newly created Fields in a named tuple;
    # Note that one major caveat to this approach is that the Fields visible to each constructor are dependent on the order
    # in which the variables were declared :/
    return foldl(vars, init=(;)) do nt, var
        # note that we call initialize here with both the current accumualated named tuple of Fields + the context given by 'fields'
        field = initialize(var, grid, clock, boundary_conditions, merge(nt, fields))
        merge(nt, (; varname(var) => field))
    end
end

"""
    initialize(var::AbstractVariable, grid::AbstractLandGrid, clock::Clock, boundary_conditions::NamedTuple, fields::NamedTuple)

Initialize a `Field` on `grid` based on the given `var` metadata. The named tuple of `boundary_conditions` should follow the standard convention of
`(var1 = (; top, bottom, ...), var2 = (; top, bottom, ...))`. If `fields` contains a `Field` matching the name of `var`, this field
will be directly returned. Otherwise, the new `Field` is constructed using the given `boundary_conditions` with the other `fields` being
made available to the constructor for auxiliary variables.
"""
function initialize(var::AbstractVariable, grid::AbstractLandGrid, clock::Clock, boundary_conditions::NamedTuple, fields::NamedTuple)
    if hasproperty(fields, varname(var))
        return getproperty(fields, varname(var))
    else
        bcs = get(boundary_conditions, varname(var), nothing)
        field = Field(grid, vardims(var), bcs)
        # if field is an input variable and has a default value/initializer, call set! on it
        if isa(var, InputVariable) && !isnothing(var.default)
            set!(field, var.default)
        end
        return field
    end
end

# Intialization for auxiliary variables that may define custom Field constructors
function initialize(var::AuxiliaryVariable, grid::AbstractLandGrid, clock::Clock, boundary_conditions::NamedTuple, fields::NamedTuple)
    if hasproperty(fields, varname(var))
        return getproperty(fields, varname(var))
    elseif isnothing(var.ctor)
        # retrieve boundary condition (if any) and create Field
        bcs = get(boundary_conditions, varname(var), nothing)
        return Field(grid, vardims(var), bcs)
    else
        # invoke field constructor if specified
        return var.ctor(grid, clock, fields)
    end
end

# Base overrides

function Adapt.adapt_structure(to, vars::StateVariables{NF}) where {NF}
    return StateVariables(
        NF,
        Adapt.adapt_structure(to, vars.prognostic),
        Adapt.adapt_structure(to, vars.tendencies),
        Adapt.adapt_structure(to, vars.auxiliary),
        Adapt.adapt_structure(to, vars.inputs),
        Adapt.adapt_structure(to, vars.namespaces),
        Adapt.adapt_structure(to, vars.clock),
    )
end

Base.eltype(::StateVariables{NF}) where {NF} = NF

Base.propertynames(
    vars::StateVariables{NF, prognames, tendnames, auxnames, inputnames, nsnames}
) where {NF, prognames, tendnames, auxnames, inputnames, nsnames} = (
    prognames...,
    auxnames...,
    inputnames...,
    nsnames...,
    fieldnames(typeof(vars))...,
)

function Base.getproperty(
    vars::StateVariables{NF, prognames, tendnames, auxnames, inputnames, nsnames},
    name::Symbol
) where {NF, prognames, tendnames, auxnames, inputnames, nsnames}
    # forward getproperty calls to variable groups
    if name ∈ prognames
        return getproperty(getfield(vars, :prognostic), name)
    elseif name ∈ auxnames
        return getproperty(getfield(vars, :auxiliary), name)
    elseif name ∈ inputnames
        return getproperty(getfield(vars, :inputs), name)
    elseif name ∈ nsnames
        return getproperty(getfield(vars, :namespaces), name)
    else
        return getfield(vars, name)
    end
end

# helper function e.g. for usage with Enzyme 
function Base.fill!(
    state::StateVariables{NF, prognames, tendnames, auxnames, namespaces}, 
    value
) where {NF, prognames, tendnames, auxnames, namespaces}
    
    for progname in prognames
        fill!(getproperty(state.prognostic, progname), value)
    end
    for tendname in tendnames
        fill!(getproperty(state.tendencies, tendname), value)
    end
    for auxname in auxnames
        fill!(getproperty(state.auxiliary, auxname), value)
    end
    return nothing 
end

function Base.copyto!(
    state::StateVariables{NF, prognames, tendnames, auxnames, inputnames, nsnames}, 
    other::StateVariables{NF, prognames, tendnames, auxnames, inputnames, nsnames}
) where {NF, prognames, tendnames, auxnames, inputnames, nsnames}
    
    for progname in prognames
        copyto!(getproperty(state.prognostic, progname), getproperty(other.prognostic, progname))
    end
    for tendname in tendnames
        copyto!(getproperty(state.tendencies, tendname), getproperty(other.tendencies, tendname))
    end
    for auxname in auxnames
        copyto!(getproperty(state.auxiliary, auxname), getproperty(other.auxiliary, auxname))
    end
    for inputname in inputnames
        copyto!(getproperty(state.inputs, inputname), getproperty(other.inputs, inputname))
    end
    for nsname in nsnames
        copyto!(getproperty(state.namespaces, nsname), getproperty(other.namespaces, nsname))
    end
    set!(state.clock, other.clock)
    return nothing
end

function Base.summary(vars::StateVariables{NF}) where {NF}
    clockstr = summary(vars.clock)
    str = "StateVariables{$NF}(clock = $clockstr, prognostic = $(keys(vars.prognostic)), auxiliary = $(keys(vars.auxiliary)), inputs = $(keys(vars.inputs)), namespaces = $(keys(vars.namespaces)))"
    return str
end

function Base.show(io::IO, vars::StateVariables{NF}) where {NF}
    println(io, "StateVariables{$NF}")
    print(io, "├─ Clock: $(vars.clock)")
    println(io)
    print(io, "├─ Prognostic: ")
    show(io, vars.prognostic)
    println(io)
    print(io, "├─ Auxiliary: ")
    show(io, vars.auxiliary)
    println(io)
    print(io, "├─ Input: ")
    show(io, vars.inputs)
    println(io)
    print(io, "├─ Namespaces: $(keys(vars.namespaces))")
end
