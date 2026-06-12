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
        prognames, closurenames, auxnames, inputnames, nsnames, cachenames,
        ProgFields, TendFields, AuxFields, InputFields, Namespaces,
        CacheFields,
        ClockType,
    } <: AbstractStateVariables
    prognostic::NamedTuple{prognames, ProgFields}
    tendencies::NamedTuple{prognames, TendFields}
    auxiliary::NamedTuple{auxnames, AuxFields}
    inputs::NamedTuple{inputnames, InputFields}
    namespaces::NamedTuple{nsnames, Namespaces}
    timestepper_cache::NamedTuple{cachenames, CacheFields}
    clock::ClockType

    function StateVariables(
            ::Type{NF},
            closurenames::Tuple{Vararg{Symbol}},
            prognostic::NamedTuple{prognames, ProgFields},
            tendencies::NamedTuple{prognames, TendFields},
            auxiliary::NamedTuple{auxnames, AuxFields},
            inputs::NamedTuple{inputnames, InputFields},
            namespaces::NamedTuple{nsnames, Namespaces},
            timestepper_cache::NamedTuple{cachenames, CacheFields},
            clock::ClockType,
        ) where {
            NF, prognames, auxnames, inputnames, nsnames, cachenames,
            ProgFields, TendFields, AuxFields, InputFields, Namespaces, CacheFields, ClockType,
        }
        return new{
            NF, prognames, closurenames, auxnames, inputnames, nsnames, cachenames,
            ProgFields, TendFields, AuxFields, InputFields, Namespaces, CacheFields, ClockType,
        }(
            prognostic,
            tendencies,
            auxiliary,
            inputs,
            namespaces,
            timestepper_cache,
            clock,
        )
    end

end

# Name getters (always type-stable, inlined constant propagation)
@inline prognostic_names(state::StateVariables) = keys(getfield(state, :prognostic))
@inline auxiliary_names(state::StateVariables) = keys(getfield(state, :auxiliary))
@inline input_names(state::StateVariables) = keys(getfield(state, :inputs))
@inline namespace_names(state::StateVariables) = keys(getfield(state, :namespaces))
@inline closure_names(::StateVariables{NF, pnames, cnames}) where {NF, pnames, cnames} = cnames

@inline timestepper_cache_names(state::StateVariables) = keys(getfield(state, :timestepper_cache))

# Allow reconstruction from properties
ConstructionBase.constructorof(::Type{StateVariables{NF, pnames, cnames}}) where {NF, pnames, cnames} = (args...) -> StateVariables(NF, cnames, args...)

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
    return if compute_tendencies
        compute_tendencies!(state, model)
    end
end

"""
Invoke `fill_halo_regions!` for all prognostic `Field`s in `state`.
"""
function Oceananigans.BoundaryConditions.fill_halo_regions!(state::StateVariables)
    # fill_halo_regions! for all prognostic variables
    fastiterate(prognostic_names(state)) do progname
        fill_halo_regions!(getproperty(state.prognostic, progname), state.clock, state)
    end

    # fill_halo_regions! for all closure variables (stored in state.auxiliary)
    fastiterate(closure_names(state)) do closurename
        fill_halo_regions!(getproperty(state.auxiliary, closurename), state.clock, state)
    end

    # recurse over namespaces
    return fastiterate(state.namespaces) do ns
        fill_halo_regions!(ns)
    end
end

"""
Reset all `Field`s in `state` to zero.
"""
function Oceananigans.TimeSteppers.reset!(state::StateVariables)
    # reset all prognostic fields
    fastiterate(state.prognostic) do field
        set!(field, zero(eltype(field)))
    end
    fastiterate(state.auxiliary) do field
        # TODO: technically we should apply auxiliary variable initializers here
        isa(field, Field) && set!(field, zero(eltype(field)))
    end
    # reset all tendency fields
    fastiterate(state.tendencies) do field
        set!(field, zero(eltype(field)))
    end
    # recurse over namespaces
    return fastiterate(state.namespaces) do ns
        reset!(ns)
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
    return fastiterate(state.namespaces) do ns
        reset_tendencies!(ns)
    end
end

"""
Initialize input variables from the given input `sources`. The `scope` corresponds to the
path of namespace names from the root namespace to `state` and is used to match namespaced
input sources to their target variables; see [`inputpath`](@ref).
"""
function initialize!(state::StateVariables, sources::InputSources, scope::Tuple{Vararg{Symbol}} = ())
    # initialize inputs in current namespace
    initialize!(state.inputs, sources, state.clock, scope)
    # recursively initialize namespaces
    fastiterate(namespace_names(state)) do nsname
        initialize!(getproperty(getfield(state, :namespaces), nsname), sources, (scope..., nsname))
    end
    return nothing
end

"""
Update input variables from the given input `sources`. The `scope` corresponds to the
path of namespace names from the root namespace to `state` and is used to match namespaced
input sources to their target variables; see [`inputpath`](@ref).
"""
function update_inputs!(state::StateVariables, sources::InputSources, scope::Tuple{Vararg{Symbol}} = ())
    # update inputs in current namespace
    update_inputs!(state.inputs, sources, state.clock, scope)
    # recursively update namespaces
    fastiterate(namespace_names(state)) do nsname
        update_inputs!(getproperty(getfield(state, :namespaces), nsname), sources, (scope..., nsname))
    end
    return
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
        else
            isa(query, Pair)
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
@inline get_field(state, var::AbstractVariable{name}) where {name} = getproperty(state, name)

"""
    $TYPEDSIGNATURES

Retrieves all `Field`s from `state` matching the names of the given variables.
"""
@inline function get_fields(state, vars::Tuple{Vararg{AbstractVariable}})
    vars = deduplicate_vars(vars)
    matched_fields = fastmap(vars) do var
        get_field(state, var)
    end
    names = map(varname, vars)
    return NamedTuple{names}(matched_fields)
end

"""
    $SIGNATURES

Retrieves all non-tendency `Field`s from `state` defined on the given `components`.
"""
@inline function get_fields(state, components...; except = (;))
    component_vars = fastmap(components) do component
        allvars = variables(component)
        closurevars = closure_variables(component)
        tuplejoin(allvars, closurevars)
    end
    vars = tuplejoin(component_vars...)
    return ntdiff(get_fields(state, vars), except)
end

"""
    $SIGNATURES

Retrieves all `Field`s from `state` corresponding to prognostic variables defined on the given `components`.
"""
@inline function prognostic_fields(state, components...)
    component_progvars = fastmap(prognostic_variables, components)
    progvars = tuplejoin(component_progvars...)
    return get_fields(state, progvars)
end

"""
    $SIGNATURES

Retrieves all `Field`s from `state` corresponding to tendencies defined on the given `components`.
"""
@inline function tendency_fields(state, components...)
    component_progvars = fastmap(prognostic_variables, components)
    progvars = tuplejoin(component_progvars...)
    return get_fields(state.tendencies, progvars)
end

"""
    $SIGNATURES

Retrieves all `Field`s from `state` corresponding to auxiliary variables defined on the given `components`.
"""
@inline function auxiliary_fields(state, components...)
    component_auxvars = fastmap(auxiliary_variables, components)
    auxvars = tuplejoin(component_auxvars...)
    return get_fields(state, auxvars)
end

"""
    $SIGNATURES

Retrieves all `Field`s from `state` corresponding to closure variables defined on the given `components`.
"""
@inline function closure_fields(state, components...)
    component_closurevars = fastmap(closure_variables, components)
    closurevars = tuplejoin(component_closurevars...)
    return get_fields(state, closurevars)
end

"""
    $SIGNATURES

Retrieves all `Field`s from `state` corresponding to input variables defined on the given `components`.
"""
@inline function input_fields(state, components...)
    component_inputvars = fastmap(input_variables, components)
    inputvars = tuplejoin(component_inputvars...)
    return get_fields(state, inputvars)
end

# Initialization of StateVariables from models and processes

"""
    $TYPEDSIGNATURES

Initialize a `StateVariables` data structure containing `Field`s for all variables defined by `model` defined on its
associated `grid`. The `clock` specifies the initial simulation time and is mutated on each time step. User-specified
`boundary_conditions` and `initializers` can be provided as `NamedTuple`s with keys corresponding to the names of state
variables to which they should be applied. If the state variables are defined within namespaces, the given `NamedTuple`
must follow the same structure. The `fields` argument allows for manual preconstruction of `Field`s for the named state
variables.
"""
function initialize(
        model::AbstractModel{NF},
        params = nothing;
        clock = Clock(time = zero(NF)),
        input_variables = (),
        timestepper = ForwardEuler{NF}(),
        boundary_conditions = (;),
        initializers = (;),
        fields = (;)
    ) where {NF}
    model_rec = isnothing(params) ? model : ParameterEditing.reconstruct(model, params)
    vars = Variables(tuplejoin(variables(model_rec), input_variables))
    state = initialize(vars, model_rec.grid; clock, timestepper, boundary_conditions, initializers, fields)
    return state
end

"""
    $TYPEDSIGNATURES

Initialize a `StateVariables` data structure containing `Field`s defined on the given `grid`
for all variables defined by `process`. Any predefined `boundary_conditions` and `fields` will
be passed through to `initialize` for each variable.
"""
function initialize(
        process::AbstractProcess{NF},
        grid::AbstractLandGrid{NF},
        params = nothing;
        clock = Clock(time = zero(NF)),
        input_variables = (),
        timestepper = ForwardEuler{NF}(),
        boundary_conditions = (;),
        initializers = (;),
        fields = (;)
    ) where {NF}
    process_rec = isnothing(params) ? process : ParameterEditing.reconstruct(process, params)
    vars = Variables(tuplejoin(variables(process_rec), input_variables))
    state = initialize(vars, grid; clock, timestepper, boundary_conditions, initializers, fields)
    return state
end

# Initialization from variable metadata

"""
    $TYPEDSIGNATURES

Initialize a `StateVariables` data structure containing `Field`s defined on the given `grid`
for all variables in `vars`. Any predefined `boundary_conditions` and `fields` will be passed
through to `initialize` for each variable.
"""
function initialize(
        @nospecialize(vars::Variables),
        grid::AbstractLandGrid{NF};
        clock::Clock = Clock(time = 0.0),
        timestepper = ForwardEuler{NF}(),
        boundary_conditions = (;),
        initializers = (;),
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
        initialize(ns.vars, grid; clock, boundary_conditions = ns_bcs, fields = ns_fields)
    end
    # get closure variable names
    closurenames = map(varname, closure_variables(values(vars.prognostic)))
    # construct StateVariables with an empty cache; the timestepper-specific cache
    # is allocated below now that all other state variables have been initialized
    initial_state = StateVariables(
        NF,
        closurenames,
        prognostic_fields,
        tendency_fields,
        auxiliary_fields,
        input_fields,
        namespaces,
        (;),
        clock,
    )
    cache = initialize(timestepper, initial_state)
    state = StateVariables(
        NF,
        closurenames,
        prognostic_fields,
        tendency_fields,
        auxiliary_fields,
        input_fields,
        namespaces,
        cache,
        clock,
    )
    # Apply Field initializers
    initialize!(state, initializers)
    return state
end

# Base case: empty named tuples
initialize(::NamedTuple{(), Tuple{}}, ::AbstractLandGrid, ::Clock, ::NamedTuple, ::NamedTuple) = (;)

"""
    $TYPEDSIGNATURES

Initialize `Field`s on `grid` for each of the variables in the given named tuple `vars`.
Any predefined `boundary_conditions` and `fields` will be passed through to `initialize`
for each variable.
"""
function initialize(
        @nospecialize(vars::NamedTuple{names, <:Tuple{Vararg{AbstractVariable}}}),
        grid::AbstractLandGrid,
        clock::Clock,
        boundary_conditions::NamedTuple,
        fields::NamedTuple
    ) where {names}
    # Initialize or retrieve Fields for each variable in `var`, accumulating the newly created Fields in a named tuple;
    # Note that one major caveat to this approach is that the Fields visible to each constructor are dependent on the order
    # in which the variables were declared :/
    return foldl(vars, init = (;)) do nt, var
        # note that we call initialize here with both the current accumualated named tuple of Fields + the context given by 'fields'
        field = initialize(var, grid, clock, boundary_conditions, merge(nt, fields))
        merge(nt, (; varname(var) => field))
    end
end

"""
    $TYPEDSIGNATURES    

Initialize a `Field` on `grid` based on the given `var` metadata. The named tuple of `boundary_conditions` should follow the standard convention of
`(var1 = (; top, bottom, ...), var2 = (; top, bottom, ...))`. If `fields` contains a `Field` matching the name of `var`, this field
will be directly returned. Otherwise, the new `Field` is constructed using the given `boundary_conditions` with the other `fields` being
made available to the constructor for auxiliary variables.
"""
function initialize(@nospecialize(var::AbstractVariable), grid::AbstractLandGrid, clock::Clock, boundary_conditions::NamedTuple, fields::NamedTuple)
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
"""
    $TYPEDSIGNATURES

Initialize a `Field` on `grid` for the given [`AuxiliaryVariable`](@ref).
"""
function initialize(@nospecialize(var::AuxiliaryVariable), grid::AbstractLandGrid, clock::Clock, boundary_conditions::NamedTuple, fields::NamedTuple)
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

function Adapt.adapt_structure(to, state::StateVariables{NF}) where {NF}
    return StateVariables(
        NF,
        closure_names(state),
        Adapt.adapt_structure(to, state.prognostic),
        Adapt.adapt_structure(to, state.tendencies),
        Adapt.adapt_structure(to, state.auxiliary),
        Adapt.adapt_structure(to, state.inputs),
        Adapt.adapt_structure(to, state.namespaces),
        Adapt.adapt_structure(to, state.timestepper_cache),
        Adapt.adapt_structure(to, state.clock),
    )
end

Base.eltype(::StateVariables{NF}) where {NF} = NF

Base.propertynames(state::StateVariables) = tuplejoin(
    prognostic_names(state),
    auxiliary_names(state),
    input_names(state),
    namespace_names(state),
    fieldnames(typeof(state)),
)

function Base.getproperty(state::StateVariables, name::Symbol)
    # forward getproperty calls to variable groups
    if name ∈ prognostic_names(state)
        return getproperty(getfield(state, :prognostic), name)
    elseif name ∈ auxiliary_names(state)
        return getproperty(getfield(state, :auxiliary), name)
    elseif name ∈ input_names(state)
        return getproperty(getfield(state, :inputs), name)
    elseif name ∈ namespace_names(state)
        return getproperty(getfield(state, :namespaces), name)
    else
        return getfield(state, name)
    end
end

# helper function e.g. for usage with Enzyme
function Base.fill!(state::StateVariables{NF}, value) where {NF}
    for progname in prognostic_names(state)
        fill!(getproperty(state.prognostic, progname), NF(value))
    end
    for tendname in keys(state.tendencies)
        fill!(getproperty(state.tendencies, tendname), NF(value))
    end
    for auxname in keys(state.auxiliary)
        fill!(getproperty(state.auxiliary, auxname), NF(value))
    end
    return nothing
end

function Base.copyto!(state::SV, other::SV) where {SV <: StateVariables}
    for progname in prognostic_names(state)
        copyto!(getproperty(state.prognostic, progname), getproperty(other.prognostic, progname))
    end
    for tendname in prognostic_names(state)
        copyto!(getproperty(state.tendencies, tendname), getproperty(other.tendencies, tendname))
    end
    for auxname in auxiliary_names(state)
        copyto!(getproperty(state.auxiliary, auxname), getproperty(other.auxiliary, auxname))
    end
    for inputname in input_names(state)
        copyto!(getproperty(state.inputs, inputname), getproperty(other.inputs, inputname))
    end
    for nsname in namespace_names(state)
        copyto!(getproperty(state.namespaces, nsname), getproperty(other.namespaces, nsname))
    end
    set!(state.clock, other.clock)
    return nothing
end

function Base.summary(state::StateVariables{NF}) where {NF}
    clockstr = summary(state.clock)
    str = "StateVariables{$NF}(clock = $clockstr, prognostic = $(keys(state.prognostic)), auxiliary = $(keys(state.auxiliary)), inputs = $(keys(state.inputs)), namespaces = $(keys(state.namespaces)), timestepper_cache = $(timestepper_cache_names(state)))"
    return str
end

function Base.show(io::IO, state::StateVariables{NF}) where {NF}
    println(io, "StateVariables{$NF}")
    print(io, "├─ Clock: $(state.clock)")
    println(io)
    print(io, "├─ Prognostic: ")
    show(io, state.prognostic)
    println(io)
    print(io, "├─ Auxiliary: ")
    show(io, state.auxiliary)
    println(io)
    print(io, "├─ Inputs: ")
    show(io, state.inputs)
    println(io)
    print(io, "├─ Namespaces: $(keys(state.namespaces))")
    println(io)
    return print(io, "└─ Timestepper cache: $(timestepper_cache_names(state))")
end
