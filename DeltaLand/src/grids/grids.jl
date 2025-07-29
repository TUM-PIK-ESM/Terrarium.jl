abstract type AbstractLandGrid end

include("vertical_discretization.jl")

# TODO: Replace with custom implementation of Oceananigans grid later
struct GlobalRingGrid{
    VertSpacing<:AbstractVerticalSpacing,
    RingGrid<:RingGrids.AbstractGrid,
    GridImpl<:Oceananigans.AbstractGrid,
} <: AbstractLandGrid
    vert::VertSpacing
    ringrid::RingGrid
    gridimpl::GridImpl
    function GlobalRingGrid(
        vert::AbstractVerticalSpacing,
        ringgrid::RingGrids.AbstractGrid
    )
        Nh = RingGrids.get_npoints(ringgrid)
        Nz = get_npoints(vert)
        # TODO: Need to consider ordering of array dimensions;
        # using the z-axis here probably results in inefficient memory access patterns
        # since most or all land computations will be along this axis
        z_coords = -cumsum(get_spacing(vert))
        gridimpl = RectilinearGrid(arch, size=(Nh, Nz), x=(0, 1), z=z_coords, topology=(Periodic, Flat, Bounded))
        return new{typeof(vert),typeof(rings),typeof(gridimpl)}(vert, ringgrid, gridimpl)
    end
end

# this is a temporary workaround to 
get_grid_impl(grid::GlobalRingGrid) = grid.gridimpl

struct OceananigansGrid{GridImpl<:Oceananigans.AbstractGrid}
    gridimpl::GridImpl
end

get_grid_impl(grid::OceananigansGrid) = grid.gridimpl
