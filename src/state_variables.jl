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
        prognames, closurenames, timestepper_classes, auxnames, inputnames, nsnames, cachenames,
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
            timestepper_classes::Tuple{Vararg{Symbol}},
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
            NF, prognames, closurenames, timestepper_classes, auxnames, inputnames, nsnames, cachenames,
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
@inline timestepper_classes(::StateVariables{NF, pnames, cnames, pclasses}) where {NF, pnames, cnames, pclasses} = pclasses

@inline timestepper_cache_names(state::StateVariables) = keys(getfield(state, :timestepper_cache))

# Allow reconstruction from properties
ConstructionBase.constructorof(::Type{StateVariables{NF, pnames, cnames, pclasses}}) where {NF, pnames, cnames, pclasses} = (args...) -> StateVariables(NF, cnames, pclasses, args...)

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
Initialize input variables from the given input `sources`.
"""
function initialize!(state::StateVariables, sources::InputSources)
    # initialize inputs in current namespace
    initialize!(state.inputs, sources, state.clock)
    # recursively initialize namespaces
    for ns in state.namespaces
        initialize!(ns, sources)
    end
    return nothing
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
# initialize model state
state = StateVariables(model)
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
# TODO: Renamed directly to `StateVariables` because without the timestepper in the `ModelIntegrator` it would be ambiguous
"""
    $TYPEDSIGNATURES

Initialize a `StateVariables` data structure containing `Field`s for all variables defined by `model` defined on its
associated `grid`. The `clock` specifies the initial simulation time and is mutated on each time step. User-specified
`boundary_conditions` and `initializers` can be provided as `NamedTuple`s with keys corresponding to the names of state
variables to which they should be applied. If the state variables are defined within namespaces, the given `NamedTuple`
must follow the same structure. The `fields` argument allows for manual preconstruction of `Field`s for the named state
variables. The time stepper cache(s) are allocated from the model's `timesteppers`. The `timestepper_classes` keyword,
a `NamedTuple` of `varname => class`, overrides the default timestepper class (`:explicit`/`:implicit`) of the named
prognostic variables (see [`prognostic`](@ref) for how defaults are declared).
"""
function StateVariables(
        model::AbstractModel{NF};
        clock = Clock(time = zero(NF)),
        input_variables = (),
        timestepper_classes = (;),
        boundary_conditions = (;),
        initializers = (;),
        fields = (;)
    ) where {NF}
    vars = Variables(tuplejoin(variables(model), input_variables))
    state = StateVariables(vars, model.grid; clock, timesteppers = get_timesteppers(model), timestepper_classes, boundary_conditions, initializers, fields)
    return state
end

"""
    $TYPEDSIGNATURES

Initialize a `StateVariables` data structure containing `Field`s defined on the given `grid`
for all variables defined by `process`. Any predefined `boundary_conditions` and `fields` will
be passed through to `initialize` for each variable.
"""
function StateVariables(
        process::AbstractProcess{NF},
        grid::AbstractLandGrid{NF};
        clock = Clock(time = zero(NF)),
        input_variables = (),
        timesteppers = default_timesteppers(NF),
        timestepper_classes = (;),
        boundary_conditions = (;),
        initializers = (;),
        fields = (;)
    ) where {NF}
    vars = Variables(tuplejoin(variables(process), input_variables))
    state = StateVariables(vars, grid; clock, timesteppers, timestepper_classes, boundary_conditions, initializers, fields)
    return state
end

# Initialization from variable metadata

"""
    $TYPEDSIGNATURES

Initialize a `StateVariables` data structure containing `Field`s defined on the given `grid`
for all variables in `vars`. Any predefined `boundary_conditions` and `fields` will be passed
through to `initialize` for each variable.

By default, each prognostic variable is integrated by the timestepper class (`:explicit` or `:implicit`)
it was declared with via [`prognostic`](@ref). The `timestepper_classes` keyword, a `NamedTuple` of
`varname => class`, overrides this default for exactly the named variables.
"""
function StateVariables(
        @nospecialize(vars::Variables),
        grid::AbstractLandGrid{NF};
        clock::Clock = Clock(time = 0.0),
        timesteppers = default_timesteppers(NF),
        timestepper_classes = (;),
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
        ns_classes = get(timestepper_classes, varname(ns), (;))
        StateVariables(ns.vars, grid; clock, timestepper_classes = ns_classes, boundary_conditions = ns_bcs, fields = ns_fields)
    end
    # get closure variable names
    closurenames = map(varname, closure_variables(values(vars.prognostic)))
    # Resolve the timestepper class (:explicit/:implicit) of each prognostic variable, in the same order as
    # the names. The default class is the one each variable was declared with (see `prognostic`); the
    # `timestepper_classes` keyword overrides it per-variable. The resolved classes are stored as a type
    # parameter so that `prognostic_names(state, class)` can select per class in a type-stable way.
    resolved_classes = consolidate_timestepper_classes(vars.prognostic, vars.namespaces, timestepper_classes)
    # construct StateVariables with an empty cache; the timestepper-specific cache
    # is allocated below now that all other state variables have been initialized
    initial_state = StateVariables(
        NF,
        closurenames,
        resolved_classes,
        prognostic_fields,
        tendency_fields,
        auxiliary_fields,
        input_fields,
        namespaces,
        (;),
        clock,
    )
    # allocate a cache for each of the model's timesteppers, keyed by timestepper class (e.g. `explicit`, `implicit`)
    cache = map(ts -> initialize(ts, initial_state), timesteppers)
    state = StateVariables(
        NF,
        closurenames,
        resolved_classes,
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

"""
    $SIGNATURES

Determine the timestepper class (`:explicit` or `:implicit`) for each prognostic variable in `progvars`,
returned as a tuple in the same order as `keys(progvars)`. Each variable's default class is the one it was
declared with via [`prognostic`](@ref); the `overrides` `NamedTuple` (mapping `varname => class`) replaces
the class for exactly the named variables. Unrecognized keys raise an `ArgumentError`.
"""
function consolidate_timestepper_classes(progvars::NamedTuple, namespaces::NamedTuple, overrides::NamedTuple)
    # validate override keys: each must name a prognostic variable at this level or a namespace
    # (namespace entries are nested NamedTuples forwarded to the namespace's own StateVariables)
    valid_keys = (keys(progvars)..., keys(namespaces)...)
    for key in keys(overrides)
        key in valid_keys || throw(
            ArgumentError(
                "`timestepper_classes` has unknown key :$key; expected a prognostic variable or namespace name in $(valid_keys)"
            )
        )
    end
    # resolve each prognostic variable's class, applying overrides and validating the result
    return map(values(progvars)) do var
        class = get(overrides, varname(var), timestepper(var)) # take the override if it exists, otherwise default to timestepper(var)
        class in (:explicit, :implicit) || throw(
            ArgumentError(
                "timestepper class for prognostic variable :$(varname(var)) must be :explicit or :implicit, got :$class"
            )
        )
        return class
    end
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
        timestepper_classes(state),
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
