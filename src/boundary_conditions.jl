abstract type AbstractBoundaryConditions end

"""
    $SIGNATURES

Constructs a suitable Oceananigans `BoundaryCondition` for the given state variable `var` on `grid`.
If `bcs` is a `NamedTuple`, it is assumed that the keys correspond to variable names and this method
is invoked recursively on the entry matching the name of `var`, if it exists. Otherwise, `nothing` is
returned.
"""
get_field_boundary_conditions(bcs::AbstractBoundaryConditions, grid::AbstractLandGrid, var::AbstractVariable) = nothing
get_field_boundary_conditions(bc::BoundaryCondition, ::AbstractLandGrid, var::AbstractVariable) = bc
function get_field_boundary_conditions(bcs::NamedTuple, grid::AbstractLandGrid, var::AbstractVariable)
    if haskey(bcs, varname(var))
        return get_field_boundary_conditions(getproperty(bcs, varname(var)), grid, var)
    else
        return nothing
    end
end

"""
Like models/processes, boundary conditions can define state variables which may be computed from
other state variables or from input data in `compute_auxiliary!`.
"""
variables(::AbstractBoundaryConditions) = ()

"""
Updates any state variables associated with the given boundary conditions.
"""
compute_auxiliary!(state, model, ::AbstractBoundaryConditions) = nothing

"""
    $SIGNATURES

Computes the boudnary tendency for the grid cell at `loc`; zero for all other grid cells.
"""
boundary_tendency(i, j, k, grid, loc, state, bc::AbstractBoundaryConditions, args...) = zero(eltype(grid))

"""
Represents default boundary conditions for all state variables. This typically
implies zero-flux boundary conditions on bounded domains and periodic boundary
conditions on periodic domains.
"""
struct DefaultBoundaryConditions <: AbstractBoundaryConditions end

"""
    $TYPEDEF

Represents a flux prescribed at the boundary of a spatial domain and applied as a forcing term.

Properties:
$TYPEDFIELDS
"""
struct PrescribedFlux{F, D<:VarDims} <: AbstractBoundaryConditions
    "Name of the flux state variable"
    name::Symbol

    "Constant, `Field`, or function corresponding to the prescribed boundary flux"
    value::F

    "Dimensions of the flux variable"
    dims::D
end

PrescribedFlux(name::Symbol, value, dims=XY()) = PrescribedFlux(name, value, dims)

variables(bc::PrescribedFlux) = (
    auxiliary(bc.name, bc.dims),
)

function compute_auxiliary!(state, model, bc::PrescribedFlux)
    F = getproperty(state, bc.name)
    set!(F, bc.value)
end

@inline function boundary_tendency(i, j, k, grid, loc, state, bc::PrescribedFlux)
    field_grid = get_field_grid(grid)
    Q = getproperty(state, bc.name)
    Δz = Δzᵃᵃᶜ(i, j, k, field_grid)
    return all(map(==, (i, j, k), loc)) * (Q[i, j, k] / Δz)
end

"""
    $TYPEDEF

Represents a the boundary conditions applied at the top and bottom of a 1D column model
discrietized along the vertical (depthwise) axis.

Properties:
$TYPEDFIELDS
"""
@kwdef struct VerticalBoundaryConditions{
    TopBC,
    BottomBC
} <: AbstractBoundaryConditions
    "Boundary condition(s) applied at the top of the vertical column."
    top::TopBC = DefaultBoundaryConditions()
    
    "Boundary condition(s) applied at the botom of the vertical column."
    bottom::BottomBC = DefaultBoundaryConditions()
end

function get_field_boundary_conditions(bcs::VerticalBoundaryConditions, grid::AbstractLandGrid, var::AbstractVariable)
    top = get_field_boundary_conditions(bcs.top, grid, var)
    bottom = get_field_boundary_conditions(bcs.bottom, grid, var)
    return field_boundary_conditions(grid; top, bottom)
end

variables(bcs::VerticalBoundaryConditions) = tuplejoin(variables(bcs.top), variables(bcs.bottom))

function compute_auxiliary!(state, model, bcs::VerticalBoundaryConditions)
    compute_auxiliary!(state, model, bcs.top)
    compute_auxiliary!(state, model, bcs.bottom)
end

"""
    field_boundary_conditions(grid::AbstractLandGrid, loc::Tuple; at...)

Creates a regularized `FieldBoundaryConditions` type from the given keyword arugments of Oceananigans `BoundaryCondition`s
with keys corresponding to their positions on the domain (i.e. `top`, `bottom`, etc.), as well as the grid and location `loc`.
The location refers the position on the staggered grid at which the boundary conditions are defined. For 1D (vertical) domains,
this is usually `(Center(), Center(), nothing)`.
"""
function field_boundary_conditions(grid::AbstractLandGrid, loc::Tuple=(Center(), Center(), nothing); at...)
    field_grid = get_field_grid(grid)
    # select only the field boundary condition arguments that are not nothing;
    # Note that since the keyword arguments are pairs we have to first select the value with `last`
    bcs = filter(!isnothing ∘ last, at)
    # create the FieldBoundaryConditions type
    field_bcs = FieldBoundaryConditions(field_grid, loc; bcs...)
    # return the regularized boundary conditions
    return regularize_field_boundary_conditions(field_bcs, field_grid, loc)
end
