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
@kwdef struct StateVariables{
    prognames, tendnames, auxnames, inputnames, nsnames, closurenames,
    ProgFields, TendFields, AuxFields, InputFields, Namespaces, Closures,
    ClockType,
} <: AbstractStateVariables
    prognostic::NamedTuple{prognames, ProgFields} = (;)
    tendencies::NamedTuple{tendnames, TendFields} = (;)
    auxiliary::NamedTuple{auxnames, AuxFields} = (;)
    inputs::NamedTuple{inputnames, InputFields} = (;)
    namespaces::NamedTuple{nsnames, Namespaces} = (;)
    closures::NamedTuple{closurenames, Closures} = (;)
    clock::ClockType = Clock()
end

function StateVariables(
    model::AbstractModel,
    clock::Clock,
    inputs::InputFields,
)
    # extract abstract variables from model
    vars = Variables(variables(model)...)
    # create named tuples for each variable type
    prognostic_vars = (; map(var -> varname(var) => var, vars.prognostic)...)
    tendency_vars = (; map(var -> varname(var) => var, vars.tendencies)...)
    auxiliary_vars = (; map(var -> varname(var) => var, vars.auxiliary)...)
    input_vars = (; map(var -> varname(var) => var, vars.inputs)...)
    namespace_vars = (; map(var -> varname(var) => var, vars.namespaces)...)
    # get grid and boundary conditions
    grid = get_grid(model)
    bcs = get_boundary_conditions(model)
    field_bcs = get_field_boundary_conditions(bcs, grid)
    # create fields from abstract variables
    init(var::AbstractVariable) = Field(grid, vardims(var), get(field_bcs, varname(var), nothing))
    # input variables are retrieved (or allocated) in the external storage provided
    init(var::InputVariable) = get_input_field(inputs, varname(var), vardims(var))
    prognostic_fields = map(init, prognostic_vars)
    tendency_fields =  map(init, tendency_vars)
    auxiliary_fields = map(init, auxiliary_vars)
    input_fields = map(init, input_vars)
    # recursively initialize state variables for each namespace
    namespaces = map(ns -> StateVariables(getproperty(model, varname(ns)), clock, inputs), namespace_vars)
    # get named tuple mapping prognostic variabels to their respective closure relations, if defined
    closures = map(var -> var.closure, filter(hasclosure, prognostic_vars))
    # construct and return StateVariables
    return StateVariables(
        prognostic_fields,
        tendency_fields,
        auxiliary_fields,
        input_fields,
        namespaces,
        closures,
        clock,
    )
end

# adapt_structure dispatch for GPU compat
function Adapt.adapt_structure(to, sv::StateVariables)
    return StateVariables(
        Adapt.adapt_structure(to, sv.prognostic),
        Adapt.adapt_structure(to, sv.tendencies),
        Adapt.adapt_structure(to, sv.auxiliary),
        Adapt.adapt_structure(to, sv.inputs),
        Adapt.adapt_structure(to, sv.namespaces),
        Adapt.adapt_structure(to, sv.closures),
        Adapt.adapt_structure(to, sv.clock),
    )
end

Base.propertynames(
    vars::StateVariables{prognames, tendnames, auxnames, inputnames, nsnames, closures}
) where {prognames, tendnames, auxnames, inputnames, nsnames, closures} = (
    prognames...,
    tendnames...,
    auxnames...,
    inputnames...,
    nsnames...,
    closures...,
    fieldnames(typeof(vars))...,
)

function Base.getproperty(
    vars::StateVariables{prognames, tendnames, auxnames, inputnames, nsnames, closures},
    name::Symbol
) where {prognames, tendnames, auxnames, inputnames, nsnames, closures}
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

function fill_halo_regions!(state::StateVariables)
    # fill_halo_regions! for all prognostic variables
    fastiterate(state.prognostic) do var
        fill_halo_regions!(var, state.clock, state)
    end
    # fill_halo_regions! for all auxiliary variables
    fastiterate(state.auxiliary) do var
        fill_halo_regions!(var, state.clock, state)
    end
    # recurse over namespaces
    fastiterate(state.namespaces) do ns
        fill_halo_regions!(ns)
    end
end

function reset_tendencies!(state::StateVariables)
    # reset all tendency fields
    for var in state.tendencies
        set!(var, zero(eltype(var)))
    end
    # recurse over namespaces
    for ns in state.namespaces
        reset_tendencies!(ns)
    end
end

# helper function e.g. for usage with Enzyme 
function Base.fill!(
    state::StateVariables{prognames, tendnames, auxnames, namespaces, closures}, 
    value
) where {prognames, tendnames, auxnames, namespaces, closures}
    
    for progname in prognames
        fill!(getproperty(state, progname), value)
    end
    for tendname in tendnames
        fill!(getproperty(state, tendname), value)
    end
    for auxname in auxnames
        fill!(getproperty(state, auxname), value)
    end
    return nothing 
end

"""
    prognostic_fields(state::StateVariables, ::Val{with_closures}=Val{true}())

Return all prognostic `Field`s in `state` in a (possibly nested) named tuple. By default,
all prognostic variables are returned (i.e. `with_closures=true`); passing `Val{false}()`
will instead return the fields for only closure-free prognostic variables.
"""
@generated function prognostic_fields(
    state::StateVariables{prognames, tendnames, auxnames, inputnames, nsnames, closures},
    ::Val{with_closures}=Val{true}()
) where {prognames, tendnames, auxnames, inputnames, nsnames, closures, with_closures}
    progn = with_closures ? prognames : filter(∉(closures), prognames)
    pv_exprs = map(pv -> :($(QuoteNode(pv)) => state.prognostic.$pv), progn)
    ns_exprs = map(ns -> :($(QuoteNode(ns)) => prognostic_fields(state.namespaces.$ns, Val{with_closures}())), nsnames)
    return :((; $(pv_exprs...), $(ns_exprs...)))
end

"""
    tendency_fields(state::StateVariables)

Return all tendency `Field`s in `state` in a (possibly nested) named tuple.
"""
@generated function tendency_fields(
    state::StateVariables{prognames, tendnames, auxnames, inputnames, nsnames}
) where {prognames, tendnames, auxnames, inputnames, nsnames}
    tv_exprs = map(tv -> :($(QuoteNode(tv)) => state.prognostic.$tv), tendnames)
    ns_exprs = map(ns -> :($(QuoteNode(ns)) => tendency_fields(state.namespaces.$ns)), nsnames)
    return :((; $(tv_exprs...), $(ns_exprs...)))
end

"""
    tendency_fields(state::StateVariables)

Return all closure `Field`s and relations in `state` in a (possibly nested) named tuple. The keys
are the prognostic variable names and the values the `Field`s of the corresponding closure variables.
"""
@generated function closure_fields(
    state::StateVariables{prognames, tendnames, auxnames, inputnames, nsnames, closures}
) where {prognames, tendnames, auxnames, inputnames, nsnames, closures}
    cv_exprs = map(pv -> :($(QuoteNode(pv)) => state.auxiliary[varname(getvar(state.closures.$pv))]), closures)
    ns_exprs = map(ns -> :($(QuoteNode(ns)) => closure_fields(state.namespaces.$ns)), nsnames)
    return :((; $(cv_exprs...), $(ns_exprs...)))
end
