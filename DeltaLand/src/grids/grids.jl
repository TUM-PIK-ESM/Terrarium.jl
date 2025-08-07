# TODO: These grid types should be replaced with proper implementations of Oceananigans `AbstractGrid`
# at some point. However, for an intiail prototype, we can just wrap `RectilinearGrid` from Oceananigans.

abstract type AbstractLandGrid{NF} end

Base.eltype(::AbstractLandGrid{NF}) where {NF} = NF

"""
    ColumnGrid{NF,RectGrid<:Grids.RectilinearGrid} <: AbstractLandGrid

Represents a set of laterally independent vertical columns with dimensions (x, y, z)
where `x` is the column dimension, `y=1` is constant, and `z` is the vertical axis.
"""
struct ColumnGrid{NF,RectGrid<:Grids.RectilinearGrid} <: AbstractLandGrid{NF}
    "Underlying Oceananigans rectilinear grid on which `Field`s are defined."
    grid::RectGrid

    """
        $SIGNATURES

    Creates a `ColumnGrid` from the given `architecture`, numeric type `NF`, and vertical
    discretization `vert`. The `num_columns` determines the number of laterally independent
    columns initialized on the grid.
    """
    function ColumnGrid(
        arch::AbstractArchitecture,
        ::Type{NF},
        vert::AbstractVerticalSpacing,
        num_columns::Int = 1,
    ) where {NF<:AbstractFloat}
        Nz = get_npoints(vert)
        # TODO: Need to eventually consider ordering of array dimensions;
        # using the z-axis here probably results in inefficient memory access patterns
        # since most or all land computations will be along this axis
        z_thick = get_spacing(vert)
        z_coords = vcat(-reverse(cumsum(z_thick)), zero(eltype(z_thick)))
        grid = Grids.RectilinearGrid(arch, size=(num_columns, Nz), x=(0, 1), z=z_coords, topology=(Periodic, Flat, Bounded))
        return new{NF,typeof(grid)}(grid)
    end
    # Default constructors
    ColumnGrid(vert::AbstractVerticalSpacing{NF}, num_columns::Int=1) where {NF} = ColumnGrid(CPU(), NF, vert, num_columns)
    ColumnGrid(arch::AbstractArchitecture, vert::AbstractVerticalSpacing{NF}, num_columns::Int=1) where {NF} = ColumnGrid(arch, NF, vert, num_columns)
    ColumnGrid(grid::Grids.RectilinearGrid{NF}) where {NF} = new{NF, typeof(grid)}(grid)
end

function Adapt.adapt_structure(to, grid::ColumnGrid)
    inner_grid = Adapt.adapt_structure(to, grid.grid)
    return ColumnGrid(inner_grid)
end

"""
    $SIGNATURES

Retrieves the underlying numerical grid on which `Field`s are defined.
"""
get_field_grid(grid::ColumnGrid) = grid.grid

"""
    $TYPEDEF

Convenience wrapper around `ColumnGrid` that defines the columns as points on a global `RingGrid`.
"""
struct GlobalRingGrid{
    NF,
    RingGrid<:RingGrids.AbstractGrid,
    RectGrid<:Grids.AbstractGrid,
} <: AbstractLandGrid{NF}
    "RingGrid specfying the lateral spatial discretization of the globe."
    rings::RingGrid

    "Underlying grid type for all independent vertical columns."
    grid::RectGrid
    
    GlobalRingGrid(vert::AbstractVerticalSpacing, rings::RingGrids.AbstractGrid) = GlobalRingGrid(CPU(), Float32, vert, rings)
    GlobalRingGrid(arch::AbstractArchitecture, vert::AbstractVerticalSpacing, rings::RingGrids.AbstractGrid) = GlobalRingGrid(arch, Float32, vert, rings)
    GlobalRingGrid(rings::RingGrid, grid::RectGrid) where {RingGrid, RectGrid} = new{eltype(grid), RingGrid, RectGrid}(rings, grid)

    """
        $SIGNATURES
    """
    function GlobalRingGrid(
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

get_field_grid(grid::GlobalRingGrid) = grid.grid

function Adapt.adapt_structure(to, grid::GlobalRingGrid)
    inner_grid = Adapt.adapt_structure(to, grid.grid)
    return GlobalRingGrid(Adapt.adapt_structure(to, grid.rings), inner_grid)
end

# Convenience dispatch for Oceananigans.launch!
function launch!(grid::AbstractLandGrid, args...; kwargs...)
    _grid = get_field_grid(grid)
    launch!(_grid.architecture, _grid, :xyz, args...; kwargs...)
end
