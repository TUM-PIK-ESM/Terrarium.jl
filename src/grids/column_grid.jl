"""
    ColumnGrid{NF,RectGrid<:OceananigansGrids.RectilinearGrid} <: AbstractLandGrid

Represents a set of laterally independent vertical columns with dimensions (x, y, z)
where `x` is the column dimension, `y=1` is constant, and `z` is the vertical axis.
"""
struct ColumnGrid{NF,RectGrid<:OceananigansGrids.RectilinearGrid} <: AbstractLandGrid{NF}
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
        grid = OceananigansGrids.RectilinearGrid(arch, size=(num_columns, Nz), x=(0, 1), z=z_coords, topology=(Periodic, Flat, Bounded))
        return new{NF,typeof(grid)}(grid)
    end
    # Default constructors
    ColumnGrid(vert::AbstractVerticalSpacing{NF}, num_columns::Int=1) where {NF} = ColumnGrid(CPU(), NF, vert, num_columns)
    ColumnGrid(arch::AbstractArchitecture, vert::AbstractVerticalSpacing{NF}, num_columns::Int=1) where {NF} = ColumnGrid(arch, NF, vert, num_columns)
    ColumnGrid(grid::OceananigansGrids.RectilinearGrid{NF}) where {NF} = new{NF, typeof(grid)}(grid)
end

get_field_grid(grid::ColumnGrid) = grid.grid

architecture(grid::ColumnGrid) = architecture(grid.grid)

function Adapt.adapt_structure(to, grid::ColumnGrid)
    inner_grid = Adapt.adapt_structure(to, grid.grid)
    return ColumnGrid(inner_grid)
end
