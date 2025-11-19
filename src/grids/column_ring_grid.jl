"""
    $TYPEDEF

Represents a global (spherical) grid of independent, vertical columns where the
spatial discretization in the horizontal direction is defined by a `RingGrids.AbstractGrid`.
"""
struct ColumnRingGrid{
        NF,
        Arch,
        RingGrid <: RingGrids.AbstractGrid,
        RectGrid <: OceananigansGrids.RectilinearGrid,
        Mask <: Union{AbstractArray, RingGrids.AbstractField},
    } <: AbstractLandGrid{NF, Arch}
    "RingGrid specfying the lateral spatial discretization of the globe"
    rings::RingGrid

    "`RingGrids.Field` (or GPU-adapted array) representing a boolean-valued mask over `rings`"
    mask::Mask

    "Underlying `Oceananigans` `RectilinearGrid` type on which `Field`s are defined"
    grid::RectGrid

    function ColumnRingGrid(
            rings::RingGrids.AbstractGrid,
            mask::AbstractArray,
            grid::OceananigansGrids.RectilinearGrid
        )
        arch = architecture(grid)
        return new{eltype(grid), typeof(arch), typeof(rings), typeof(grid), typeof(mask)}(rings, mask, grid)
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
            mask::RingGrids.AbstractField{Bool} = convert.(Bool, ones(rings))
        ) where {NF <: AbstractFloat}
        assert_field_matches_grid(mask, rings)
        # get number of horizontal grid poitns by summing over mask
        Nh = sum(mask)
        # get number of vertical grid points
        Nz = num_layers(vert)
        # TODO: Need to consider ordering of array dimensions;
        # using the z-axis here probably results in inefficient memory access patterns
        # since most or all land computations will be along this axis
        z_thick = get_spacing(vert)
        z_coords = convert.(NF, vcat(-reverse(cumsum(z_thick)), zero(eltype(z_thick))))
        grid = OceananigansGrids.RectilinearGrid(arch, NF, size = (Nh, Nz), x = (1, Nh), z = z_coords, topology = (Periodic, Flat, Bounded))
        # adapt ring grid and mask
        rings = on_architecture(arch, rings)
        mask = on_architecture(arch, mask)
        return new{NF, typeof(arch), typeof(rings), typeof(grid), typeof(mask)}(rings, mask, grid)
    end

    ColumnRingGrid(
        arch::AbstractArchitecture,
        vert::AbstractVerticalSpacing,
        rings::RingGrids.AbstractGrid,
        mask::RingGrids.AbstractField{Bool} = convert.(Bool, ones(rings))
    ) = ColumnRingGrid(arch, Float32, vert, rings, mask)

    ColumnRingGrid(
        arch::AbstractArchitecture,
        ::Type{NF},
        vert::AbstractVerticalSpacing,
        mask::RingGrids.AbstractField{Bool} = convert.(Bool, ones(rings))
    ) where {NF} = ColumnRingGrid(arch, NF, vert, mask.grid, mask)

    ColumnRingGrid(
        arch::AbstractArchitecture,
        vert::AbstractVerticalSpacing,
        mask::RingGrids.AbstractField{Bool}
    ) = ColumnRingGrid(arch, Float32, vert, mask.grid, mask)

    ColumnRingGrid(
        vert::AbstractVerticalSpacing,
        rings::RingGrids.AbstractGrid,
        mask::RingGrids.AbstractField{Bool} = convert.(Bool, ones(rings))
    ) = ColumnRingGrid(CPU(), Float32, vert, rings, mask)

    ColumnRingGrid(
        vert::AbstractVerticalSpacing,
        mask::RingGrids.AbstractField{Bool}
    ) = ColumnRingGrid(CPU(), Float32, vert, mask.grid, mask)
end

@adapt_structure ColumnRingGrid

get_field_grid(grid::ColumnRingGrid) = grid.grid

"""
    $SIGNATURES

Converts the given Oceananigans `Field` to a `RingGrids.Field` with a ring grid matching that of the given `ColumnRingGrid`.
"""
RingGrids.Field(field::Field{LX, LY, Nothing}, grid::ColumnRingGrid; fill_value = NaN) where {LX, LY} = RingGrids.Field(architecture(field), interior(field), grid; fill_value)
RingGrids.Field(field::Field{LX, LY, LZ}, grid::ColumnRingGrid; fill_value = NaN) where {LX, LY, LZ} = RingGrids.Field(architecture(field), dropdims(interior(field), dims = 2), grid; fill_value)
RingGrids.Field(field::AbstractVecOrMat, grid::ColumnRingGrid; fill_value = NaN) = RingGrids.Field(architecture(grid), field, grid; fill_value)
function RingGrids.Field(arch::AbstractArchitecture, field::AbstractVecOrMat, grid::ColumnRingGrid; fill_value = NaN)
    # need to be on CPU to do the copying
    grid = on_architecture(arch, grid)
    field = on_architecture(arch, field)
    # create new RingGrids field initialized with fill_value
    ring_field = RingGrids.Field(grid.rings, size(field)[2:end]...)
    fill!(ring_field, fill_value)
    # need to access underlying data arrays directly to avoid scalar indexing
    colons = (Colon() for d in size(field)[2:end])
    ring_field.data[grid.mask.data, colons...] .= field
    return ring_field
end

function on_architecture(arch::AbstractArchitecture, grid::ColumnRingGrid)
    return ColumnRingGrid(
        on_architecture(arch, grid.rings),
        on_architecture(arch, grid.mask),
        on_architecture(arch, grid.grid)
    )
end

function Base.show(io::IO, mime::MIME"text/plain", grid::ColumnRingGrid{NF}) where {NF}
    println(io, "ColumnRingGrid{$NF} on $(architecture(grid)) with")
    show(io, mime, grid.rings)
    println(io)
    return show(io, mime, grid.grid)
end
