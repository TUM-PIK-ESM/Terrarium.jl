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

    function StateVariables(
        ::Type{NF},
        prognostic::NamedTuple{prognames, ProgFields},
        tendencies::NamedTuple{tendnames, TendFields},
        auxiliary::NamedTuple{auxnames, AuxFields},
        inputs::NamedTuple{inputnames, InputFields},
        namespaces::NamedTuple{nsnames, Namespaces},
        closures::NamedTuple{closurenames, Closures},
        clock::ClockType
    ) where {NF, prognames, tendnames, auxnames, inputnames, nsnames, closurenames,
             ProgFields, TendFields, AuxFields, InputFields, Namespaces, Closures, ClockType}
        return new{NF, prognames, tendnames, auxnames, inputnames, nsnames, closurenames,
                   ProgFields, TendFields, AuxFields, InputFields, Namespaces, Closures, ClockType}(
            prognostic,
            tendencies,
            auxiliary,
            inputs,
            namespaces,
            closures,
            clock
        )
    end
end

# Allow reconstruction from properties
ConstructionBase.constructorof(::Type{StateVariables{NF}}) where {NF} = (args...) -> StateVariables(NF, args...)

function StateVariables(
    model::AbstractModel{NF},
    clock::Clock,
    inputs::InputFields,
) where {NF}
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
        NF,
        prognostic_fields,
        tendency_fields,
        auxiliary_fields,
        input_fields,
        namespaces,
        closures,
        clock,
    )
end

function Adapt.adapt_structure(to, vars::StateVariables{NF}) where {NF}
    return StateVariables(
        NF,
        Adapt.adapt_structure(to, vars.prognostic),
        Adapt.adapt_structure(to, vars.tendencies),
        Adapt.adapt_structure(to, vars.auxiliary),
        Adapt.adapt_structure(to, vars.inputs),
        Adapt.adapt_structure(to, vars.namespaces),
        Adapt.adapt_structure(to, vars.closures),
        Adapt.adapt_structure(to, vars.clock),
    )
end

Base.eltype(::StateVariables{NF}) where {NF} = NF

Base.propertynames(
    vars::StateVariables{NF, prognames, tendnames, auxnames, inputnames, nsnames, closures}
) where {NF, prognames, tendnames, auxnames, inputnames, nsnames, closures} = (
    prognames...,
    auxnames...,
    inputnames...,
    nsnames...,
    closures...,
    fieldnames(typeof(vars))...,
)

function Base.getproperty(
    vars::StateVariables{NF, prognames, tendnames, auxnames, inputnames, nsnames, closures},
    name::Symbol
) where {NF, prognames, tendnames, auxnames, inputnames, nsnames, closures}
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
    state::StateVariables{NF, prognames, tendnames, auxnames, namespaces, closures}, 
    value
) where {NF, prognames, tendnames, auxnames, namespaces, closures}
    
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
    state::StateVariables{NF, prognames, tendnames, auxnames, inputnames, nsnames, closurenames}, 
    other::StateVariables{NF, prognames, tendnames, auxnames, inputnames, nsnames, closurenames}
) where {NF, prognames, tendnames, auxnames, inputnames, nsnames, closurenames}
    
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
    # currently clock isn't copied, not defined, and not our type
    #copyto!(state.clock, other.clock)
    return nothing
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
Call `invclosure!` (i.e. inverse mapping from closure to prognostic variables) on
all closure relations defined in `state`.
"""
function invclosure!(state::StateVariables, model::AbstractModel)
    for closure in state.closures
        invclosure!(state, model, closure)
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
