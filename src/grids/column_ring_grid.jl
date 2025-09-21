"""
    $TYPEDEF

Convenience wrapper around `ColumnGrid` that defines the columns as points on a global `RingGrid`.
"""
struct ColumnRingGrid{
    NF,
    RingGrid<:RingGrids.AbstractGrid,
    RectGrid<:Grids.AbstractGrid,
} <: AbstractLandGrid{NF}
    "RingGrid specfying the lateral spatial discretization of the globe."
    rings::RingGrid

    "Underlying grid type for all independent vertical columns."
    grid::RectGrid
    
    ColumnRingGrid(vert::AbstractVerticalSpacing, rings::RingGrids.AbstractGrid) = ColumnRingGrid(CPU(), Float32, vert, rings)
    ColumnRingGrid(arch::AbstractArchitecture, vert::AbstractVerticalSpacing, rings::RingGrids.AbstractGrid) = ColumnRingGrid(arch, Float32, vert, rings)
    ColumnRingGrid(rings::RingGrid, grid::RectGrid) where {RingGrid, RectGrid} = new{eltype(grid), RingGrid, RectGrid}(rings, grid)

    """
        $SIGNATURES
    """
    function ColumnRingGrid(
        arch::AbstractArchitecture,
        ::Type{NF},
        vert::AbstractVerticalSpacing,
        rings::RingGrids.AbstractGrid
    ) where {NF<:AbstractFloat}
        Nh = RingGrids.get_npoints(rings)
        Nz = get_npoints(vert)
        # TODO: Need to consider ordering of array dimensions;
        # using the z-axis here probably results in inefficient memory access patterns
        # since most or all land computations will be along this axis
        z_thick = get_spacing(vert)
        z_coords = vcat(-reverse(cumsum(z_thick)), zero(eltype(z_thick)))
        grid = Grids.RectilinearGrid(arch, NF, size=(Nh, Nz), x=(0, 1), z=z_coords, topology=(Periodic, Flat, Bounded))
        return new{NF,typeof(rings),typeof(grid)}(rings, grid)
    end
end

get_field_grid(grid::ColumnRingGrid) = grid.grid

function Adapt.adapt_structure(to, grid::ColumnRingGrid)
    inner_grid = Adapt.adapt_structure(to, grid.grid)
    return ColumnRingGrid(Adapt.adapt_structure(to, grid.rings), inner_grid)
end
