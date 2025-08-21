abstract type AbstractBoundaryConditions end

"""
Alias for a `NamedTuple` of `BoundaryCondition` types.
"""
const VarBoundaryConditions{names, BCs} = NamedTuple{names, BCs} where {names, BCs<:Tuple{Vararg{BoundaryCondition}}}

"""
    $SIGNATURES

Constructs a suitable Oceananigans `BoundaryCondition` for the given state variable `var` on `grid`.
If `bcs` is a `NamedTuple`, it is assumed that the keys correspond to variable names and this method
is invoked recursively on the entry matching the name of `var`, if it exists. Otherwise, `nothing` is
returned.
"""
get_field_boundary_conditions(::AbstractBoundaryConditions, grid::AbstractLandGrid) = (;)
get_field_boundary_conditions(bc::BoundaryCondition, ::AbstractLandGrid) = bc
get_field_boundary_conditions(bcs::VarBoundaryConditions, ::AbstractLandGrid) = bcs

"""
Like models/processes, boundary conditions can define state variables which may be computed from
other state variables or from input data in `compute_auxiliary!`.
"""
variables(::AbstractBoundaryConditions) = ()
variables(bc::BoundaryCondition) = ()
variables(bcs::VarBoundaryConditions) = tuplejoin(map(variables, bcs)...)

"""
Updates any state variables associated with the given boundary conditions.
"""
compute_auxiliary!(state, model, ::AbstractBoundaryConditions) = nothing
compute_auxiliary!(state, model, ::BoundaryCondition) = nothing
compute_auxiliary!(state, model, ::VarBoundaryConditions) = nothing

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
struct PrescribedFlux{name, F, D<:VarDims} <: AbstractBoundaryConditions
    "Constant, `Field`, or function corresponding to the prescribed boundary flux"
    value::F

    "Dimensions of the flux variable"
    dims::D

    PrescribedFlux(name::Symbol, value, dims=XY()) = new{name, typeof(value), typeof(dims)}(value, dims)
end

@inline varname(::PrescribedFlux{name}) where {name} = name

variables(bc::PrescribedFlux) = (
    auxiliary(varname(bc), bc.dims),
)

function compute_auxiliary!(state, model, bc::PrescribedFlux)
    F = getproperty(state, varname(bc))
    set!(F, bc.value)
end

@inline function boundary_tendency(i, j, k, grid, loc, state, bc::PrescribedFlux)
    field_grid = get_field_grid(grid)
    Q = getproperty(state, varname(bc))
    Δz = Δzᵃᵃᶜ(i, j, k, field_grid)
    return all(map(==, (i, j, k), loc)) * (Q[i, j, k] / Δz)
end

"""
    $TYPEDEF

Represents the boundary conditions applied at the top and bottom of a 1D column model
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

function get_field_boundary_conditions(
    bcs::VerticalBoundaryConditions,
    grid::AbstractLandGrid,
)
    top = get_field_boundary_conditions(bcs.top, grid)
    bottom = get_field_boundary_conditions(bcs.bottom, grid)
    # invert the nested structure of the top/bottom named tuples;
    # this way we have a single named tuple of variable names where the entries 
    var_bcs = merge_recursive(map(bc -> (top=bc,), top), map(bc -> (bottom=bc,), bottom))
    return map(var_bcs) do bc
        FieldBoundaryConditions(grid, (Center(), Center(), nothing); bc...)
    end
end

variables(bcs::VerticalBoundaryConditions) = tuplejoin(variables(bcs.top), variables(bcs.bottom))

function compute_auxiliary!(state, model, bcs::VerticalBoundaryConditions)
    compute_auxiliary!(state, model, bcs.top)
    compute_auxiliary!(state, model, bcs.bottom)
end

"""
    FieldBoundaryConditions(grid::AbstractLandGrid, loc::Tuple; at...)

Creates a regularized `FieldBoundaryConditions` type from the given keyword arugments of Oceananigans `BoundaryCondition`s
with keys corresponding to their positions on the domain (i.e. `top`, `bottom`, etc.), as well as the grid and location `loc`.
The location refers the position on the staggered grid at which the boundary conditions are defined. For 1D (vertical) domains,
this is usually `(Center(), Center(), nothing)`.
"""
function FieldBoundaryConditions(grid::AbstractLandGrid, loc::Tuple; at...)
    field_grid = get_field_grid(grid)
    bcs = map(((k, bc)) -> k => isnothing(bc) ? NoFluxBoundaryCondition() : bc, keys(at), values(at))
    # create the FieldBoundaryConditions type
    field_bcs = FieldBoundaryConditions(field_grid, loc; immersed=DefaultBoundaryCondition(), bcs...)
    # return the regularized boundary conditions
    return regularize_field_boundary_conditions(field_bcs, field_grid, loc)
end
