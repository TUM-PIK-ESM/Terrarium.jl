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
    ClockType,
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
    namespace_vars = filter(var -> isa(var, Namespace), vars)
    # get tendencies from prognostic variables
    tendencies = map(var -> var.tendency, prognostic)
    # create closure variables and add to auxiliary variables
    closurevars = map(var -> getvar(var.closure, vardims(var)), filter(hasclosure, prognostic))
    auxiliary_ext = tuplejoin(auxiliary, closurevars)
    # merge all variable names
    varnames = tuplejoin(map(varname, prognostic), map(varname, auxiliary_ext), map(varname, namespace_vars))
    @assert merge_duplicates(varnames) == varnames "all state variable names within one namespace must be unique! found duplicates in $(varnames)"
    # get progvar => closure pairs
    closures = map(var -> varname(var) => var.closure, filter(hasclosure, prognostic))
    # get grid, initializer, and boundary conditions
    grid = get_grid(model)
    init = get_initializer(model)
    bcs = get_boundary_conditions(model)
    # intialize fields
    prognostic_fields = map(var -> varname(var) => create_field(var, init, bcs, grid), prognostic)
    tendency_fields = map(var -> varname(var) => create_field(var, init, bcs, grid), tendencies)
    auxiliary_fields = map(var -> varname(var) => create_field(var, init, bcs, grid), auxiliary_ext)
    # recursively initialize state variables for namespaces
    namespaces = map(ns -> varname(ns) => StateVariables(getproperty(model, varname(ns)), clock), namespace_vars)
    # construct and return StateVariables struct
    return StateVariables(
        (; prognostic_fields...),
        (; tendency_fields...),
        (; auxiliary_fields...),
        (; namespaces...),
        (; closures...),
        clock,
    )
end

# adapt_structure dispatch for GPU compat
function Adapt.adapt_structure(to, sv::StateVariables)
    return StateVariables(
        Adapt.adapt_structure(to, sv.prognostic),
        Adapt.adapt_structure(to, sv.tendencies),
        Adapt.adapt_structure(to, sv.auxiliary),
        Adapt.adapt_structure(to, sv.namespaces),
        Adapt.adapt_structure(to, sv.closures),
        Adapt.adapt_structure(to, sv.clock),
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

function fill_halo_regions!(state::StateVariables)
    # fill_halo_regions! for all prognostic variables
    for var in state.prognostic
        fill_halo_regions!(var)
    end
    # fill_halo_regions! for all auxiliary variables
    for var in state.auxiliary
        fill_halo_regions!(var)
    end
    # recurse over namespaces
    for ns in state.namespaces
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
    # Specify BCs if defined
    boundary_conditions = get_field_boundary_conditions(bcs, grid, var)
    field = if !isnothing(boundary_conditions)
        FT(get_field_grid(grid), args...; boundary_conditions, kwargs...)
    else
        FT(get_field_grid(grid), args...; kwargs...)
    end
    # Apply initializer if defined
    field_init = get_field_initializer(init, grid, var)
    if !isnothing(field_init)
        set!(field, field_init)
    end
    return field
end
