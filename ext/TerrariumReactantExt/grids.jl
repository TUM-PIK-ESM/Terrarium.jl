# Grid transfers to/from the device.
#
# The generic land-grid `on_architecture` (src/grids/grids.jl) uses `adapt(array_type(arch), …)`,
# which for a Reactant array type strips the `Field`/grid architecture to `nothing`. Instead we
# transfer the *inner* Oceananigans grid with Oceananigans' own `on_architecture` (which preserves
# `arch = ReactantState` and produces concrete device arrays) and rewrap it.

Terrarium.on_architecture(arch::RARCH, grid::ColumnGrid) =
    ColumnGrid(on_architecture(arch, grid.grid))

Terrarium.on_architecture(arch::CPU, grid::ColumnGrid{<:Any, <:RARCH}) =
    ColumnGrid(on_architecture(arch, grid.grid))

# ColumnRingGrid: transfer the inner Oceananigans grid to the device; keep `rings`/`mask` on the
# CPU, where they are used only for host-side pre/post-processing (masked-index field conversion).
Terrarium.on_architecture(arch::RARCH, grid::ColumnRingGrid) =
    ColumnRingGrid(grid.rings, grid.mask, on_architecture(arch, grid.grid))

Terrarium.on_architecture(arch::CPU, grid::ColumnRingGrid{<:Any, <:RARCH}) =
    ColumnRingGrid(
    on_architecture(arch, grid.rings), on_architecture(arch, grid.mask),
    on_architecture(arch, grid.grid)
)

# Ring→Oceananigans field conversion for device grids. The generic method
# (src/grids/column_ring_grid.jl) transfers the ring data to the device and gathers the masked
# points there with boolean indexing, which silently returns wrong results on Reactant's concrete
# arrays. Instead, gather the masked points on the host and transfer the finished interior with
# `set!` (which detours through the CPU for device fields).
function Oceananigans.Field(
        ring_field::RingGrids.AbstractField, grid::ColumnRingGrid{<:Any, <:RARCH};
        default_value = zero(eltype(ring_field))
    )
    mask = Array(grid.mask.data)
    data = Array(ring_field.data)
    if ndims(ring_field) == 1
        # 2D field (horizontal only)
        field = Field{Center, Center, Nothing}(grid.grid)
        values = reshape(data[mask], :, 1, 1)
    elseif ndims(ring_field) == 2
        # 3D field (horizontal + vertical)
        @assert size(grid.grid, 3) == size(ring_field, 2) "Vertical dimension mismatch: grid has $(size(grid.grid, 3)) layers, but field has $(size(ring_field, 2)) layers"
        field = Field{Center, Center, Center}(grid.grid)
        values = reshape(data[mask, :], :, 1, size(ring_field, 2))
    else
        error("Unsupported number of dimensions for RingGrids.Field: $(ndims(ring_field))")
    end
    fill!(field, default_value)
    set!(field, values)
    return field
end
