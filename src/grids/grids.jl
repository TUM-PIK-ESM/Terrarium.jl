# TODO: These grid types should be replaced with proper implementations of Oceananigans `AbstractGrid`
# at some point. However, for an intiail prototype, we can just wrap `RectilinearGrid` from Oceananigans.

abstract type AbstractLandGrid{NF, Arch} end

Base.eltype(::AbstractLandGrid{NF}) where {NF} = NF
Base.summary(grid::AbstractLandGrid{NF, Arch}) where {NF, Arch} = "$(nameof(typeof(grid))){$NF, $Arch} with dimensions $(size(grid))"
Base.size(grid::AbstractLandGrid) = size(get_field_grid(grid))

"""
Return the number of vertical layers defined by the given `grid`.
"""
num_layers(grid::AbstractLandGrid) = size(get_field_grid(grid), 3)

"""
    get_field_grid(grid::AbstractLandGrid)::Oceananigans.AbstractGrid

Returns the underlying `Oceananigans` grid type for `Field`s defined on the given land `grid`.
"""
function get_field_grid end

architecture(grid::AbstractLandGrid) = architecture(get_field_grid(grid))
on_architecture(arch, grid::AbstractLandGrid) = adapt(array_type(arch), grid)

include("grid_utils.jl")

include("column_grid.jl")

include("column_ring_grid.jl")
