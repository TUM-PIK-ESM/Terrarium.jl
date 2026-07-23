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
