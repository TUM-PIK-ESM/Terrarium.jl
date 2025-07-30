abstract type AbstractLandGrid{NF} end

# TODO: Replace with custom implementation of Oceananigans grid later
struct GlobalRingGrid{
    NF,
    RingGrid<:RingGrids.AbstractGrid,
    GridImpl<:Oceananigans.AbstractGrid,
} <: AbstractLandGrid{NF}
    "RingGrid representing the lateral spatial discretization of the globe."
    rings::RingGrid

    "Oceananigans grid implementation used for finite volume discretization of the vertical axis."
    gridimpl::GridImpl
    
    GlobalRingGrid(vert::AbstractVerticalSpacing, rings::RingGrids.AbstractGrid) = GlobalRingGrid(CPU(), vert, rings)
    function GlobalRingGrid(
        arch::AbstractArchitecture,
        vert::AbstractVerticalSpacing,
        rings::RingGrids.AbstractGrid
    )
        Nh = RingGrids.get_npoints(rings)
        Nz = get_npoints(vert)
        # TODO: Need to consider ordering of array dimensions;
        # using the z-axis here probably results in inefficient memory access patterns
        # since most or all land computations will be along this axis
        z_thick = get_spacing(vert)
        z_coords = vcat(-reverse(cumsum(z_thick)), zero(eltype(z_thick)))
        gridimpl = RectilinearGrid(arch, size=(Nh, Nz), x=(0, 1), z=z_coords, topology=(Periodic, Flat, Bounded))
        return new{eltype(gridimpl),typeof(rings),typeof(gridimpl)}(rings, gridimpl)
    end
end

# this is a temporary workaround to 
get_grid_impl(grid::GlobalRingGrid) = grid.gridimpl

struct OceananigansGrid{GridImpl<:Oceananigans.AbstractGrid}
    gridimpl::GridImpl
end

get_grid_impl(grid::OceananigansGrid) = grid.gridimpl

# Convenience dispatch for Oceananigans.launch!
function launch!(grid::AbstractLandGrid, args...; kwargs...)
    grid_impl = get_grid_impl(grid)
    launch!(grid_impl.architecture, grid_impl, :xyz, args...; kwargs...)
end