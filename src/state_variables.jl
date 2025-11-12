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
    prognames, tendnames, auxnames, inputnames, nsnames, closurenames,
    ProgFields, TendFields, AuxFields, InputFields, Namespaces, Closures,
    ClockType
} <: AbstractStateVariables
    prognostic::NamedTuple{prognames, ProgFields}
    tendencies::NamedTuple{tendnames, TendFields}
    auxiliary::NamedTuple{auxnames, AuxFields}
    inputs::NamedTuple{inputnames, InputFields}
    namespaces::NamedTuple{nsnames, Namespaces}
    closures::NamedTuple{closurenames, Closures}
    clock::ClockType
end

@adapt_structure StateVariables

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
    # get closure fields
    closure_fields = map(closure -> auxiliary_fields[varname(closurevar(closure))], values(closures))
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

# TODO: just take the first prognostic variable for the eltype, bad idea?
Base.eltype(vars::StateVariables) = eltype(vars.prognostic[1])

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

function Base.copyto!(
    state::StateVariables{prognames, tendnames, auxnames, inputnames, nsnames, closurenames}, 
    other::StateVariables{prognames, tendnames, auxnames, inputnames, nsnames, closurenames}
) where {prognames, tendnames, auxnames, inputnames, nsnames, closurenames}
    
    for progname in prognames
        copyto!(getproperty(state, progname), getproperty(other, progname))
    end
    for tendname in tendnames
        copyto!(getproperty(state, tendname), getproperty(other, tendname))
    end
    for auxname in auxnames
        copyto!(getproperty(state, auxname), getproperty(other, auxname))
    end
    for inputname in inputnames
        copyto!(getproperty(state.inputs, inputname), getproperty(other.inputs, inputname))
    end
    for nsname in nsnames
        copyto!(getproperty(state.namespaces, nsname), getproperty(other.namespaces, nsname))
    end
    # currently clock isn't copied, not defined, and not our type
    #copyto!(state.clock, other.clock)
    return nothing
end

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
