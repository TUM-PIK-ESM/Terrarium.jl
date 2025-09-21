"""
    $TYPEDEF

Represents a global (spherical) grid of independent, vertical columns where the
spatial discretization in the horizontal direction is defined by a `RingGrids.AbstractGrid`.
The resulting 
"""
struct ColumnRingGrid{
    NF,
    RingGrid<:RingGrids.AbstractGrid,
    MaskField<:RingGrids.AbstractField{Bool},
    RectGrid<:OceananigansGrids.RectilinearGrid,
} <: AbstractLandGrid{NF}
    "RingGrid specfying the lateral spatial discretization of the globe"
    rings::RingGrid

    "`RingGrids.Field` representing a boolean-valued mask over `rings`"
    mask::MaskField

    "Underlying `Oceananigans` `RectilinearGrid` type on which `Field`s are defined"
    grid::RectGrid

    function ColumnRingGrid(
        rings::RingGrids.AbstractGrid,
        mask::RingGrids.AbstractField{Bool},
        grid::OceananigansGrids.RectilinearGrid
    )
        assert_field_matches_grid(mask, rings)
        new{eltype(grid), typeof(rings), typeof(mask), typeof(grid)}(rings, mask, grid)
    end

    """
        $SIGNATURES

    Constructs a `ColumnRingGrid` over the given `rings`.
    """
    function ColumnRingGrid(
        arch::AbstractArchitecture,
        ::Type{NF},
        vert::AbstractVerticalSpacing,
        rings::RingGrids.AbstractGrid,
        mask::RingGrids.AbstractField{Bool}=convert.(Bool, ones(rings))
    ) where {NF<:AbstractFloat}
        assert_field_matches_grid(mask, rings)
        # get number of horizontal grid poitns by summing over mask
        Nh = sum(mask)
        # get number of vertical grid points
        Nz = get_npoints(vert)
        # TODO: Need to consider ordering of array dimensions;
        # using the z-axis here probably results in inefficient memory access patterns
        # since most or all land computations will be along this axis
        z_thick = get_spacing(vert)
        z_coords = vcat(-reverse(cumsum(z_thick)), zero(eltype(z_thick)))
        grid = OceananigansGrids.RectilinearGrid(arch, NF, size=(Nh, Nz), x=(0, 1), z=z_coords, topology=(Periodic, Flat, Bounded))
        return new{NF, typeof(rings), typeof(mask), typeof(grid)}(rings, mask, grid)
    end

    ColumnRingGrid(
        arch::AbstractArchitecture,
        vert::AbstractVerticalSpacing,
        rings::RingGrids.AbstractGrid,
        mask::RingGrids.AbstractField{Bool}=convert.(Bool, ones(rings))
    ) = ColumnRingGrid(arch, Float32, vert, rings, mask)

    ColumnRingGrid(
        vert::AbstractVerticalSpacing,
        rings::RingGrids.AbstractGrid,
        mask::RingGrids.AbstractField{Bool}=convert.(Bool, ones(rings))
    ) = ColumnRingGrid(CPU(), Float32, vert, rings, mask)
end

get_field_grid(grid::ColumnRingGrid) = grid.grid

function Adapt.adapt_structure(to, grid::ColumnRingGrid)
    inner_grid = Adapt.adapt_structure(to, grid.grid)
    return ColumnRingGrid(Adapt.adapt_structure(to, grid.rings), inner_grid)
end
