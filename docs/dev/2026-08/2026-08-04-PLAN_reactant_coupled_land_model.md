# Reactant support for the coupled `LandModel` (soil + snow, no vegetation)

> Status: **in progress**. All six Reactant blockers found on the coupled soil+snow path are fixed and
> the configuration compiles and runs; the CPU-vs-Reactant correctness comparison is being verified.

Date of initial draft: 2026-08-04

Base revision: e015a29db9d97090edfd112a5eea165640134404

## Originating prompt

> Add a LandModel with Soil (including Richards eq) and Snow but without Vegetation to the Reactant
> model setup and test (and work on) it running correctly with Reactant. Report issues to me.

## Revision log

- 2026-08-04: initial draft, written after the blockers were found and fixed.

## Problem description

`test/reactant/` only covered *standalone* models: `SoilModel` with the default `NoFlow` hydrology
(`:soil_heat_column`, `:soil_heat_column_stretched`, `:soil_heat_global`) and a standalone `SnowModel`
(`:snow_column`). The coupled `LandModel` — soil hydrology via the Richardson-Richards equation, a
single-layer snowpack, the surface energy balance with its implicit skin-temperature solve, and the
surface hydrology that couples them — had never been traced. It did not compile.

## Background

Six independent blockers sit between a `LandModel(grid; soil, snow, vegetation = nothing)` and a
working Reactant compile. They are listed in the order they surface, since each one hides the next.

### 1. In-kernel `@assert` with an interpolated message

`saturation_to_pressure!` (`src/processes/soil/hydrology/soil_hydraulic_closures.jl`) asserted
`isfinite(ψm)` with a string-interpolated message. String formatting inside a kernel is not merely a
trap: it drags `IOBuffer`, `jl_string_to_genericmemory`, `ijl_rethrow` and a GC frame into the device
IR. All 25 `InvalidIRError` reasons in the first failure traced back to this single line.

This is the AGENTS.md "no reachable throw paths in kernels" rule. The `@assert` was also redundant: a
`debughook!` for `saturation_to_pressure_kernel!` already runs `checkfinite!` on the output, host-side,
under `TERRARIUM_DEBUG=true`.

**Fix**: drop the `@assert`, and clamp the inverse-SWRC argument from above as well as below
(`clamp(sat * por, sat_res, por)` instead of `max(sat * por, sat_res)`). Oversaturated states would
otherwise evaluate the inverse van Genuchten curve at a negative base and produce the NaN the assertion
was watching for. For legal states (`sat ≤ 1`, guaranteed by `adjust_saturation_profile!`) the clamp is
a no-op.

### 2. Complex arithmetic in the van Genuchten hydraulic conductivity

`hydraulic_conductivity(::AbstractSoilHydraulics{NF, <:VanGenuchten, UnsatKVanGenuchten}, …)` used
`complex()` to tolerate illegal state values, on the reasoning that complex `^` and `sqrt` are total
where their real counterparts raise `DomainError`. Complex math lowers to libdevice symbols that
Reactant cannot resolve:

```
JIT session error: Symbols not found: [ ___nv_cospif, ___nv_hypotf, ___nv_sinpif, ___nv_sincosf ]
CompilationError: MLIR pass pipeline "all" failed   (gpu_compute_hydraulics_kernel!)
```

**Fix**: clamp the effective saturation to the unit interval and use real arithmetic. Every fractional
power then sees a non-negative base, so the `DomainError` path the `complex` trick was avoiding is
unreachable by construction. Verified bit-identical to the complex form for legal states.

### 3. `set!(::Field, ::Field)` across a location mismatch

`initialize!(state, grid, ::ImplicitSkinTemperature, …)` seeded `skin_temperature` (a
`Center, Center, Nothing` surface field) from `ground_temperature`, which is a *window*
(`Center, Center, Center`) of the soil temperature field. `set!` dispatches on that mismatch to
Oceananigans' interpolating fallback, which detours through the CPU and calls `pointer` on a
`ConcretePJRTArray`:

```
ERROR: conversion to pointer not defined for ConcretePJRTArray{Float32, 3, 1}
```

**Fix**: copy the interiors directly. They are the same size, so a plain broadcast is both correct and
traceable.

### 4. Solvers writing a `Field` with a 2D index

Oceananigans defines `setindex!` on a `Field` for *exactly three* indices. A 2D write falls through to
`Base`'s generic `AbstractArray` path (`_to_subscript_indices` → `axes(::AbstractField)` →
`size(::AbstractGrid, loc, indices)`), which is a dynamic dispatch on the device:

```
InvalidIRError: compiling ... gpu_solve_surface_energy_balance_kernel!
Reason: unsupported call to an unknown function (call to jl_f_throw_methoderror)
```

Both `RootSolver` and `FixedPointSolver` index their target field with the objective's index tuple,
which is `(i, j)` for the `XY` skin-temperature solve. This is AGENTS.md pitfall #10; note that 2D
*reads* happen to work, because `Field` has a catch-all `getindex(f, inds...)` that forwards to the
underlying array — there is no matching `setindex!`.

**Fix**: a `pad_indices(field, indices)` helper in `src/utils/kernel_utils.jl` expands the index tuple
to `ndims(field)` with trailing ones; both solvers use it for every target-field access. It compiles to
nothing — the generated IR is instruction-for-instruction identical to a hand-written `field[i, j, 1]`.

### 5. A lazy `AbstractOperation` as a boundary condition

`LandModel` built the water boundary condition as `InfiltrationFlux(-infiltration)`. `-field` does not
compute anything; it builds a lazy `UnaryOperation`. Boundary conditions are evaluated inside kernels
via `getbc`, which indexes the condition with two indices:

```julia
@inline getbc(bc::BC{<:Any, <:AbstractArray}, i::Integer, j::Integer, args...) = @inbounds bc.condition[i, j]
```

As in #4, that works for a `Field` (catch-all `getindex`) but not for an `AbstractOperation`, which
only defines the 3-index form:

```
InvalidIRError: compiling ... gpu__compute_z_bcs!
Reason: unsupported call to an unknown function (call to jl_f_throw_methoderror)
 [2] axes  @ Oceananigans/src/Fields/abstract_field.jl:62
 [7] getbc @ Oceananigans/src/BoundaryConditions/boundary_condition.jl:182
```

**Fix**: materialize the negation into a `soil_water_flux` coupling field (upwards-positive), declared
by `interface_variables(::LandModel)` alongside `soil_heat_flux` and refreshed by
`compute_soil_water_flux!` at the end of each `compute_auxiliary!`. This makes the water BC structurally
identical to the heat BC, which already compiles. It is semantically equivalent because
`fill_halo_regions!` (where `getbc` runs) precedes `compute_auxiliary!` in `update_state!`, so the lazy
operation also saw the previous step's infiltration.

### 6. Data-dependent control flow in kernels

Two remaining constructs could not be raised to StableHLO:

- **Masked stores.** `compute_hydraulics!` (`soil_hydrology_rre.jl`) used `if k <= 1 / elseif k >= Nz /
  else` with a store in each branch, plus a second store to face `k + 1` in one branch only:
  `error: masked affine.store is dependent on less dimensions than masked stored value`. Rewritten to
  select the *indices* and issue two unconditional stores (AGENTS.md: `ifelse`, never `if`/`else`, in
  kernels).
- **Convergence-tested loops.** The default skin-temperature solver is `RootSolver`, wrapping
  RootSolvers.jl's Newton iteration. Its `while` loop has a data-dependent trip count and lowers to an
  `scf.while` that the raise pass rejects: `error: cannot raise op to stablehlo`.

**Fix**: the masked stores are fixed in place. For the solver, a new `NewtonSolver` performs a *fixed*
number of iterations carried in the type, so the loop unrolls into straight-line code. The Reactant test
configuration selects it explicitly; **the default is unchanged** (see Known limitations).

## Summary of changes

| File | Change |
| --- | --- |
| `src/processes/soil/hydrology/soil_hydraulic_closures.jl` | Drop in-kernel `@assert`; clamp the inverse-SWRC argument to `[sat_res, por]` |
| `src/processes/soil/hydrology/soil_hydraulic_properties.jl` | Replace `complex` arithmetic in the van Genuchten conductivity with a clamped real form |
| `src/processes/soil/hydrology/soil_hydrology_rre.jl` | `compute_hydraulics!`: branchless index selection, unconditional stores |
| `src/processes/surface/skin_temperature.jl` | `initialize!`: copy interiors instead of `set!` across a location mismatch |
| `src/utils/kernel_utils.jl` | New `pad_indices` helper |
| `src/solvers/root_solvers.jl`, `src/solvers/fixed_point.jl` | Index the target field in 3D via `pad_indices` |
| `src/solvers/newton.jl` (new) | `NewtonSolver`: fixed, type-level iteration count |
| `src/solvers/solvers.jl` | Include the new solver |
| `src/models/coupled/land_model.jl` | `soil_water_flux` coupling variable + `compute_soil_water_flux!`; BC no longer a lazy operation |
| `test/reactant/setup.jl`, `test/reactant/runtests.jl` | New `:land_soil_snow` configuration |
| `test/solvers/test_solvers.jl` | `NewtonSolver` tests |
| `test/coupled_models/land_model_tests.jl` | Refresh `soil_water_flux` before asserting on the BC condition |

## Testing and verification

- `test/reactant/runtests.jl` gains `:land_soil_snow`: a single column on an `ExponentialSpacing` grid,
  Richards hydrology with a van Genuchten SWRC and `UnsatKVanGenuchten` conductivity, `SingleLayerSnow`,
  no vegetation, under sub-freezing prescribed air temperature. CPU and Reactant are stepped 100 times
  and every prognostic and auxiliary field is compared.
- The pack is held below freezing deliberately: the comparison should not hinge on the exact step at
  which a melt threshold is crossed, which CPU and XLA need not agree on.
- The full CPU suite passes unchanged.

## Known limitations

- **The default `LandModel` still does not compile under Reactant.** `ImplicitSkinTemperature` defaults
  to `RootSolver`, whose convergence-tested loop is blocker #6. The Reactant configuration overrides it
  with `NewtonSolver`; users of the default configuration hit the raise failure. Making `NewtonSolver`
  the default is recommended (see Future work) but changes results for existing users and so is left for
  review.
- **`getbc`'s 2D indexing is an upstream issue.** Terrarium now avoids it by keeping BC conditions as
  `Field`s, but any operation-valued boundary condition will hit it again. Worth reporting to
  Oceananigans.
- **Vegetation is still unsupported** — see the separate root-fraction `exp(::TracedRArray)` blocker.

## Future work

- Make `NewtonSolver` the default skin-temperature solver, after checking the effect on the SEB
  regression tests.
- Add a `LandModel` configuration to `test/reactant/autodiff.jl` once the forward path is verified.
- Report the `getbc` 2D-indexing limitation upstream.
