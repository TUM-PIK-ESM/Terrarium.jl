# Boundary conditions interface

abstract type AbstractBoundaryConditions end

"""
    get_field_boundary_conditions(bcs::AbstractBoundaryConditions, var::AbstractVariable)

Retrieve the `Field` boundary conditions for the corresponding variable. Defaults to returning `nothing` if no BC is defined.
"""
get_field_boundary_conditions(::AbstractBoundaryConditions, ::AbstractVariable, ::AbstractLandGrid) = nothing

struct PrescribedFluxes <: AbstractBoundaryConditions end

struct FieldBoundaryConditions{BCS<:NamedTuple} <: AbstractBoundaryConditions
    var_bcs::BCS
end

FieldBoundaryConditions(; vars...) = FieldBoundaryConditions((; vars...))

function get_field_boundary_conditions(bcs::FieldBoundaryConditions, var::AbstractVariable, grid::AbstractLandGrid)
    bcs = get(bcs.var_bcs, varname(var), nothing)
    if isnothing(bcs)
        return nothing
    else
        field_bc = BoundaryConditions.FieldBoundaryConditions(get_field_grid(grid), (Center,Center,Nothing); bcs...)
        return _regularize_field_boundary_conditions(field_bc, get_field_grid(grid))
    end
end

# TODO: Temporary workaround to issue with Oceananigans.jl BC regularization
# Default implementation of regularize_field_boundary_conditions in Oceananigans.jl assigns an immersed boundary
# which results in a warning for our 1D RectilinearGrid.
function _regularize_field_boundary_conditions(bcs::BoundaryConditions.FieldBoundaryConditions, grid)
    loc = (Center, Center, Nothing)
    west = BoundaryConditions.regularize_west_boundary_condition(bcs.west, grid, loc, 1, BoundaryConditions.LeftBoundary,  nothing)
    east = BoundaryConditions.regularize_east_boundary_condition(bcs.east, grid, loc, 1, BoundaryConditions.RightBoundary, nothing)
    south = BoundaryConditions.regularize_south_boundary_condition(bcs.south, grid, loc, 2, BoundaryConditions.LeftBoundary,  nothing)
    north = BoundaryConditions.regularize_north_boundary_condition(bcs.north, grid, loc, 2, BoundaryConditions.RightBoundary, nothing)
    bottom = BoundaryConditions.regularize_bottom_boundary_condition(bcs.bottom, grid, loc, 3, BoundaryConditions.LeftBoundary,  nothing)
    top = BoundaryConditions.regularize_top_boundary_condition(bcs.top, grid, loc, 3, BoundaryConditions.RightBoundary, nothing)
    return BoundaryConditions.FieldBoundaryConditions(; west, east, south, north, bottom, top)
end
