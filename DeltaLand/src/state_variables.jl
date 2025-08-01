abstract type AbstractStateVariables end

# temporary solution for holding state variables
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
    init::AbstractInitializer,
    bcs::AbstractBoundaryConditions,
    grid::AbstractLandGrid,
    clock::Clock,
    vars...
)
    # check for duplicates
    varnames = map(varname, vars)
    @assert merge_duplicates(varnames) == varnames "all state variable names within one namespace must be unique! found one or more duplicates in $(varnames)"
    # filter out variables from tuple by type
    prognostic = merge_duplicates(filter(var -> isa(var, PrognosticVariable), vars))
    auxiliary = merge_duplicates(filter(var -> isa(var, AuxiliaryVariable), vars))
    namespaces = filter(var -> isa(var, Pair{Symbol}), vars)
    # get tendencies from prognostic variables
    tendencies = map(var -> var.tendency, prognostic)
    # get progvar => closure pairs
    closures = map(var -> varname(var) => var.closure, filter(hasclosure, prognostic))
    # intialize fields
    prognostic_fields = map(var -> varname(var) => create_field(var, init, bcs, grid), prognostic)
    tendency_fields = map(var -> varname(var) => create_field(var, init, bcs, grid), tendencies)
    auxiliary_fields = map(var -> varname(var) => create_field(var, init, bcs, grid), auxiliary)
    # recursively initialize state variables for namespaces
    namespaces = map(kv -> first(kv) => StateVariables(init, bcs, grid, clock, last(kv)...), namespaces)
    return StateVariables(
        prognostic=(; prognostic_fields...),
        tendencies=(; tendency_fields...),
        auxiliary=(; auxiliary_fields...),
        namespaces=(; namespaces...),
        closures=(; closures...),
        clock,
    )
end

function StateVariables(model::AbstractModel)   
    init = get_initializer(model)
    bcs = get_boundary_conditions(model)
    grid = get_grid(model)
    vars = variables(model)
    return StateVariables(init, bcs, grid, clock, vars...)
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
    for var in getfield(state, :tendencies)
        set!(var, zero(eltype(var)))
    end
end

# Field construction

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
    inits = get_field_boundary_conditions(init, var)
    if !isnothing(inits)
        Oceananigans.set!(field, inits[name])
    end
    return field
end
