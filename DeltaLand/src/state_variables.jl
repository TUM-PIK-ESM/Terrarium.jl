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
    prognames, tendnames, auxnames, subnames, closurenames,
    ProgVars, TendVars, AuxVars, SubVars, Closures,
    ClockType <: Clock,
} <: AbstractStateVariables
    prognostic::NamedTuple{prognames, ProgVars} = (;)
    tendencies::NamedTuple{tendnames, TendVars} = (;)
    auxiliary::NamedTuple{auxnames, AuxVars} = (;)
    namespaces::NamedTuple{subnames, SubVars} = (;)
    closures::NamedTuple{closurenames, Closures} = (;)
    clock::ClockType = Clock()
end

function StateVariables(
    model::AbstractModel,
    clock::Clock
)
    vars = variables(model)
    # filter out variables from tuple by type
    prognostic = merge_duplicates(filter(var -> isa(var, PrognosticVariable), vars))
    auxiliary = merge_duplicates(filter(var -> isa(var, AuxiliaryVariable), vars))
    namespaces = filter(var -> isa(var, Namespace), vars)
    # get tendencies from prognostic variables
    tendencies = map(var -> var.tendency, prognostic)
    # create closure variables and add to auxiliary variables
    closurevars = map(var -> AuxiliaryVariable(varname(var.closure), vardims(var)), filter(hasclosure, prognostic))
    auxiliary_ext = tuplejoin(auxiliary, closurevars)
    # merge all variable names
    varnames = tuplejoin(map(varname, prognostic), map(varname, auxiliary_ext), map(varname, namespaces))
    @assert merge_duplicates(varnames) == varnames "all state variable names within one namespace must be unique! found duplicates in $(varnames)"
    # get progvar => closure pairs
    closures = map(var -> varname(var) => var.closure, filter(hasclosure, prognostic))
    # get grid, boundary conditions and initializer
    grid = get_grid(model)
    init = get_initializer(model)
    bcs = get_boundary_conditions(model)
    # intialize fields
    prognostic_fields = map(var -> varname(var) => create_field(var, init, bcs, grid), prognostic)
    tendency_fields = map(var -> varname(var) => create_field(var, init, bcs, grid), tendencies)
    auxiliary_fields = map(var -> varname(var) => create_field(var, init, bcs, grid), auxiliary_ext)
    # recursively initialize state variables for namespaces
    namespace_states = map(ns -> varname(ns) => StateVariables(getproperty(model, varname(ns)), clock), namespaces)
    return StateVariables(
        prognostic=(; prognostic_fields...),
        tendencies=(; tendency_fields...),
        auxiliary=(; auxiliary_fields...),
        namespaces=(; namespace_states...),
        closures=(; closures...),
        clock=clock,
    )
end

Base.propertynames(
    vars::StateVariables{prognames, tendnames, auxnames, namespaces, closures}
) where {prognames, tendnames, auxnames, namespaces, closures} = (
    prognames...,
    tendnames...,
    auxnames...,
    namespaces...,
    closures...,
    fieldnames(typeof(vars))...,
)

function Base.getproperty(
    vars::StateVariables{prognames, tendnames, auxnames, namespaces, closures},
    name::Symbol
) where {prognames, tendnames, auxnames, namespaces, closures}
    # forward getproperty calls to variable groups
    if name ∈ prognames
        return getproperty(getfield(vars, :prognostic), name)
    elseif name ∈ tendnames
        return getproperty(getfield(vars, :tendencies), name)
    elseif name ∈ auxnames
        return getproperty(getfield(vars, :auxiliary), name)
    elseif name ∈ namespaces
        return getproperty(getfield(vars, :namespaces), name)
    else
        return getfield(vars, name)
    end
end

function reset_tendencies!(state::StateVariables)
    for var in state.tendencies
        set!(var, zero(eltype(var)))
    end
    for substate in state.namespaces
        reset_tendencies!(substate)
    end
end

# Field construction

## Retrieves the Oceananigans `Field` type for the given variable dimensions.
## The first three type parameters refer to the location of the variable on the staggered finite volume grid:
## Center refers to the grid cell centroid, Face to the boundary, and Nothing to a quantity that has no dimensionality
## in that direction.
field_type(::XY) = Field{Center,Center,Nothing}
field_type(::XYZ) = Field{Center,Center,Center}

"""
    $SIGNATURES

Initializes a `Field` on `grid` for the variable, `var`, using the given initializer, `init`, and boundary conditions, `bcs`.
Additional arguments are passed direclty to the `Field` constructor. The location of the `Field` is determined by the specified
`VarDims` on `var`.
"""
function create_field(
    var::AbstractVariable,
    init::AbstractInitializer,
    bcs::AbstractBoundaryConditions,
    grid::AbstractLandGrid,
    args...;
    kwargs...
)
    FT = field_type(vardims(var))
    bcs = get_field_boundary_conditions(bcs, var)
    field = if !isnothing(bcs)
        FT(get_field_grid(grid), args...; boundary_conditions=bcs, kwargs...)
    else
        FT(get_field_grid(grid), args...; kwargs...)
    end
    # Apply initializer if defined
    field_init = get_field_initializer(init, var)
    if !isnothing(field_init)
        set!(field, field_init)
    end
    return field
end
