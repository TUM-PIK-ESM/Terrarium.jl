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

function StateVariables(
    vars::Variables,
    grid::AbstractLandGrid{NF},
    clock::Clock;
    boundary_conditions = (;),
    fields = (;)
) where {NF}
    # Initialize Fields for each variable group, if they are not already given in the user defined `fields`
    input_fields = initialize(vars.inputs, grid, clock, boundary_conditions, fields)
    tendency_fields = initialize(vars.tendencies, grid, clock, boundary_conditions, fields)
    auxiliary_fields = initialize(vars.auxiliary, grid, clock, boundary_conditions, merge(fields, input_fields))
    prognostic_fields = initialize(vars.prognostic, grid, clock, boundary_conditions, merge(fields, input_fields, auxiliary_fields))
    # recursively initialize state variables for each namespace
    namespaces = map(vars.namespaces) do ns
        StateVariables(ns.vars, grid, clock; boundary_conditions=get(boundary_conditions, varname(ns), (;)), fields=get(fields, varname(ns), (;)))
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

"""
    update_state!(state::StateVariables, model::AbstractModel, inputs::InputSources; compute_tendencies = true)

Update the `state` for the given `model` and `inputs`; this includes calling `update_inputs!` and
`fill_halo_regions!` followed by `compute_auxiliary!` and `compute_tendencies!`, if `compute_tendencies = true`.
"""
function update_state!(state::StateVariables, model::AbstractModel, inputs::InputSources; compute_tendencies = true)
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
function fill_halo_regions!(state::StateVariables)
    # fill_halo_regions! for all prognostic variables
    fastiterate(state.prognostic) do field
        fill_halo_regions!(field, state.clock, state)
    end
    # fill_halo_regions! for all auxiliary variables
    fastiterate(state.auxiliary) do field
        fill_halo_regions!(field, state.clock, state)
    end
    # recurse over namespaces
    fastiterate(state.namespaces) do ns
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
    get_fields(state::StateVariables, queries::Union{Symbol, Pair}...)

Retrieves fields with names given in `queries` and returns them in a `NamedTuple`. Each argument
in `queries` can either be a `Symbol` corresponding to a field/variable defined in the namespace
of `state` or a `Pair{Symbol, Tuple}` where the key is the child namespace and the value is a
tuple of queries from that namespace.

```julia
# initialize model
state = initialize(model)
# get the temperature and saturation_water_ice fields
fields = get_fields(state, :temperature, :saturation_water_ice)
# extract temperature as well as variables from a namespace
nested_fields = get_fields(state, :temperature, :namespace => (:subvar1, :subvar2))
```
"""
function get_fields(state::StateVariables, queries::Union{Symbol, Pair}...)
    fields = map(queries) do query
        if isa(query, Symbol)
            query => getproperty(state, query)
        else isa(query, Pair)
            key, value = query
            key => get_fields(getproperty(state, key), value...)
        end
    end
    return (; fields...)
end

# Default initialize dispatch for model types

function initialize(
    model::AbstractModel{NF};
    clock = Clock(time=zero(NF)),
    external_variables = (),
    boundary_conditions = (;),
    fields = (;)
) where {NF}
    vars = Variables(tuplejoin(variables(model), external_variables))
    grid = get_grid(model)
    state = StateVariables(vars, grid, clock; boundary_conditions, fields)
    return state
end

function initialize(
    vars::NamedTuple{names, <:Tuple{Vararg{AbstractVariable}}},
    grid::AbstractLandGrid,
    clock::Clock,
    boundary_conditions::NamedTuple,
    fields::NamedTuple
)
    # initialize or retrieve Fields for each variable in `var`, accumulating the newly created Fields in a named tuple
    return foldl(vars, init=(;)) do nt, var
        # note that we call initialize here with both the current accumualated named tuple of Fields + the context given by 'fields'
        field = initialize(var, grid, boundary_conditions, merge(nt, fields))
        merge(nt, (; varname(var) => field))
    end
end

function initialize(var::AbstractVariable, grid::AbstractLandGrid, clock::Clock, boundary_conditions::NamedTuple, fields::NamedTuple)
    if hasproperty(fields, varname(var))
        return getproperty(fields, varname(var))
    else
        bcs = get(boundary_conditions, varname(var), nothing)
        field = Field(grid, vardims(var), bcs)
        # if field is an input variable and has an init value/function, call set! on the specified initial value
        isa(field, InputVariable) && !isnothing(field.init) && set!(field, field.init)
    end
end

function initialize(var::AuxiliaryVariable, grid::AbstractLandGrid, clock::Clock, boundary_conditions::NamedTuple, fields::NamedTuple)
    if hasproperty(fields, varname(var))
        return getproperty(fields, varname(var))
    elseif isnothing(var.ctor)
        # retrieve boundary condition (if any) and create Field
        bcs = get(boundary_conditions, varname(var), nothing)
        return Field(grid, vardims(var), bcs)
    else
        # invoke field constructor if specified
        return Field(var.ctor(grid, clock, fields))
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