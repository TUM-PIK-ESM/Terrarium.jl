# Phase 2: SimpleNonlinearSolve algorithms inside a KernelAbstractions kernel for per-column implicit timestepping

> Status: **completed**. In-kernel per-column backward-Euler solves with SimpleNonlinearSolve algorithms,
> following the established `src/solvers/root_solvers.jl` residual-handling pattern (an `ObjectiveFunction`
> wrapped by a `solve!` that owns the iterate↔field writes and invokes the root-finder). Residual assembled
> from Terrarium's scalar kernel functions (`energy_to_temperature!` + `compute_energy_tendency`), validated
> to reproduce `compute_f!` exactly both host-side and in-kernel. Example runs end-to-end: host-vs-kernel
> agreement 0.0, cross-column agreement 0.0, per-algorithm residuals |G(u*)|∞ ≤ 2e-3, kernel-launch median
> time ~5 ms for 8 columns × N = 100.

Date of initial draft: 2026-08-05

Base revision: db48d370

## Originating prompt

> Start phase 2, where you focus on the use of simplenonlinearsolve algorithms and try to make them work with the kernel context.

> Look at how in @src/solvers/root_solvers.jl the residual function is handled. This can serve as inspiration as to how to do something similar for the solve here.

## Revision log

- 2026-08-05: Initial draft after phase-1 completion.

## Problem description

Phase 1 (see `docs/dev/2026-08/2026-08-05_PLAN_nonlinearsolve_implicit_timestep.md`) established host-side
backward-Euler implicit timestepping with NonlinearSolve.jl and documented stage 2 as design notes only.
Phase 2 implements stage 2: **one SimpleNonlinearSolve solve per soil column, launched as a
KernelAbstractions kernel** (one work-item per column), with the residual handled the same way
`src/solvers/root_solvers.jl` handles scalar root solves.

## Background

### The root_solvers.jl pattern (template)

`src/solvers/root_solvers.jl` defines the canonical Terrarium way to run a per-grid-point nonlinear solve
inside a kernel:

- An `ObjectiveFunction{target, F, DF}` (`src/solvers/solvers.jl:15`) wraps a per-point residual
  `func(out, indices..., grid, fields, func_args...)` that reads the current iterate from the target field
  and returns the residual `F(x)`.
- `solve!(out, indices, grid, fields, objective, solver, func_args...)` (`root_solvers.jl:31`) is a
  `@propagate_inbounds` kernel function that:
  - calls `build_residual(...)` which returns the closure `residual!(x) { target_field[indices...] = x; step_func!(...) }`
    (i.e., it **writes the iterate `x` back into the field**, then evaluates the objective), plus an optional
    derivative closure,
  - invokes the actual root-finder (`RootSolvers.find_zero`) on that residual,
  - writes the converged root back into the target field.
- This `solve!` is already invoked from within a `@kernel` in production
  (`surface_energy_balance.jl:141` → `solve_skin_temperature!` → `solve!`).

### Validation performed (probes)

- **All 6 SimpleNonlinearSolve algorithms** (`SimpleKlement`, `SimpleDFSane`, `SimpleLimitedMemoryBroyden`,
  `SimpleBroyden`, `SimpleNewtonRaphson`, `SimpleTrustRegion`) compile, run, and are **allocation-free**
  (0 bytes host-side with StaticArrays) and give correct results inside a KA CPU `@kernel` on a 6-vector
  test problem (dev vs host ≤ 4e-14).
- **Per-column residual assembly**: the terrain scalar kernel functions
  `energy_to_temperature!` (pointwise closure) and `compute_energy_tendency` (conduction stencil) called
  per-k over a column's `fields` NamedTuple reproduce `compute_f!` **exactly**:
  - host-side manual scalar loop: max |Δ| = 0.0 (at init and at a perturbed state exercising other freeze-curve branches);
  - **in-kernel** (one work-item per column, 3 columns): max |Δ| = 0.0 vs `compute_f!`.
- The residual needs no tendency buffer: with NoFlow hydrology and carbon-only biogeochemistry,
  `compute_energy_tendency` is the sole contributor (`compute_tendencies_kernel!` only does `tend += compute_energy_tendency`),
  so per-k `G(u) = u − u_prev − Δt·f_k(u)` with `f_k = compute_energy_tendency(...)` reproduces the residual exactly.

### Kernel-relevant facts established in phase 1 / probes

- Frozen-ghost semantics: `electrical_temperature` ghost cells are frozen at initialization
  (top = 2·T_bc − T(N_z), rest zero); `compute_f!` never re-fills halos, and the in-kernel residual
  reproduces this by only writing interior cells while leaving the Field's halo cells untouched.
- `@index` usage: the `@kernel` macro rewrites a bare `lhs = @index(...)` (defaults to Linear, returns the
  linear Int for a 1-D column dimension); parenthesized or `NamedModule.@index` forms are **not** rewritten.
- SimpleNonlinearSolve is bundled inside NonlinearSolve as `NonlinearSolve.SimpleNonlinearSolve`.
  `SciMLBase.ImmutableNonlinearProblem` is reached as `NonlinearSolve.SciMLBase.ImmutableNonlinearProblem`.
- `SimpleKlement.__solve` accepts `abstol`, `reltol`, `maxiters` and is allocation-free on StaticArrays.

## Summary of changes

1. New example `examples/autodiff/implicit_timestepping_simplenonlinearsolve_kernel.jl`
   (`julia --project=examples ...`), a Literate-style experiment that:
   - sets up a `ColumnGrid(CPU(), FT, UniformSpacing(), N_COLS)` with `N_COLS` laterally independent columns;
   - extracts the `fields`/`out` NamedTuples and physics args exactly as the soil energy closure/
     tendency path does (`Terrarium.get_fields`, `closure_fields`);
   - defines a **vector** `ObjectiveFunction` (`Terrarium.ObjectiveFunction`) wrapping
     `backward_euler_column_residual!(out, i, j, grid, fields, u_prev, Δt, ::Val{N}, func_args...)`
     assembled from `energy_to_temperature!` + `compute_energy_tendency` per k;
   - defines `solve_implicit_column!` mirroring `root_solvers.jl`'s `solve!`: builds the
     `residual!(x)` closure (writes x into `fields.internal_energy`, calls the objective), constructs an
     `ImmutableNonlinearProblem{false}`, calls `solve(prob, alg; abstol, reltol, maxiters)` with a
     SimpleNonlinearSolve algorithm, and writes the converged `sol.u` back into the field;
   - launches a `@kernel` with one work-item per column.
   - verifies each in-kernel column solution against a host SimpleNonlinearSolve reference (same residual)
     and the Ariadne `newton_krylov!` reference; benchmarks.
2. No `src/` changes (experiment-only); the design is written to mirror `src/solvers/` so it can be lifted
   into production after sign-off.

## Testing and verification

- In-kernel per-column residual == `compute_f!` (max |Δ| = 0.0).
- In-kernel column solutions agree with the host SimpleNonlinearSolve solution for the same residual.
- In-kernel column solutions agree with the Ariadne reference (≤ ~2.5e-2 absolute, as in phase 1).
- All successful algorithms reported; per-algorithm iteration counts/residual metadata where available.

## Documentation changes

- The example file itself (Literate comments) documents usage, the residual-handling pattern, and findings.
- This plan doc tracks status.

## Known limitations

- CPU only in this phase; GPU/Reactant (kernel args must become `isbits` data buffers rather than
  Oceananigans `Field` structs) is follow-up.
- StaticArrays `SVector{N_z}` with N_z = 100 per work-item raises register pressure; fine on CPU,
  to be re-evaluated on GPU (possibly smaller N_z or a block-based decomposition).
- SimpleNonlinearSolve's per-column solve uses StaticArrays `\`/LU internally only in
  `SimpleNewtonRaphson`/`SimpleTrustRegion` (unreliable for N ≳ 14); the factorization-free algorithms
  (`SimpleKlement`, `SimpleDFSane`, `SimpleLimitedMemoryBroyden`) are the primary targets.

## Future work

- Lift the pattern into `src/solvers/` as a production solver type (e.g. `ImplicitTimestepSolver`) once
  sign-off and GPU compatibility are confirmed.
- GPU route: pass `fields` as `isbits` data buffers (plain arrays over a column, indices fixed at (i,j))
  instead of Oceananigans `Field` structs; verify Reactant.
- Differentiable in-kernel residual (Enzyme/ForwardDiff over the scalar kernel functions) for AD through
  the implicit solve.
- Generalize to multiple tendency contributors (`tend += f`), e.g. hydrology + carbon, once those are added.

---

# Phase 3: Host-side Newton solver with ClimaTimeSteppers + Enzyme Jacobian

> Status: **completed**. Host-side backward-Euler Newton solver using ClimaTimeSteppers.jl `NewtonsMethod`
> for the nonlinear solve and Enzyme forward-mode JVPs (column-by-column via `Duplicated(state, dstate)`)
> for the Jacobian. The residual is computed through Fields (not flat vectors) so Enzyme can trace the
> full U→T→tendencies→residual chain. Validated: forward vs backward Euler agree to ~1.7% relative error
> after 10 timesteps with Δt = 300s.

Date of initial draft: 2026-08-06

Base revision: (current)

## Originating prompt

> Wire CTS Newton as a Terrarium timestepper. Use CTS only for the nonlinear solve, not for the ODE
> interface. Use Enzyme for the Jacobian. Store residual and Jacobian as Fields (like tendencies) for
> Enzyme compatibility.

## Summary of changes

### New file: `src/timesteppers/backward_euler.jl`

- `BackwardEuler{NF}` struct: wraps `ClimaTimeSteppers.NewtonsMethod` + `Δt`
- `timestepping(::BackwardEuler) = Implicit()` — integrates with the trait dispatch in `model_integrator.jl`
- `compute_f!(state, grid, soil, constants)`: combined `fill_halo_regions!` + `reset_tendencies!` +
  `closure!` + `compute_tendencies!` (frozen-ghost aware)
- `backward_euler_residual!(res, u, u_prev, state, grid, soil, constants, Δt)`: computes CTS-compatible
  residual `f_CTS(x) = x_prev + Δt·f(x) − x` through Fields, then extracts to flat vector for CTS
- `compute_jacobian!(J, u, u_prev, dstate, state, grid, soil, constants, Δt)`: computes Jacobian column-by-column
  using `Enzyme.autodiff(Forward, compute_f!, Const, Duplicated(state, dstate), Const(grid), ...)`
  following the pattern from `examples/autodiff/differentiating_terrarium.jl`
- `timestep!` method: orchestrates the solve — copies state to flat vector, allocates shadow state,
  builds closures, calls `CTS.solve_newton!`, writes solution back

### Key design decisions

1. **Residual through Fields**: The residual computation goes through Fields (not just flat vectors) so
   Enzyme can trace the U→T→tendencies→residual chain. Flat-vector conversion only at the CTS interface.
2. **Jacobian via column-by-column JVPs**: `Enzyme.jacobian(Forward, ...)` fails with
   `EnzymeMutabilityException` because the closure captures mutable state. The working pattern is
   `Enzyme.autodiff(Forward, compute_f!, Const, Duplicated(state, dstate), ...)` with one-hot seeding
   per column — same pattern used in `differentiating_terrarium.jl`.
3. **CTS used minimally**: Only `NewtonsMethod`, `allocate_cache`, `solve_newton!` from ClimaTimeSteppers.
   No `ODEFunction`, `ClimaODEFunction`, or `ODEProblem`.
4. **`vec(interior(field))` not `parent(field)`**: Excludes halos; uses Oceananigans-recommended access pattern.
   `copy()` is required because CTS mutates `x` in-place via `x .-= Δx`.
5. **`I` (UniformScaling) cannot be broadcast**: Use explicit `J[i,i] -= one(eltype(u))` loop instead of
   `J .= Δt .* J_f .- I`.

### Dependencies added to `Project.toml`

- `ClimaTimeSteppers = "4ffaf3f8-f178-474f-9d79-c298d935aac7"` (v0.10.5)
- `Enzyme = "7da31134-29cb-4562-af13-18bfc0fc98df"` (v0.13.198)
- `LinearAlgebra` (stdlib)

### Imports added to `src/Terrarium.jl`

- `import ClimaTimeSteppers`
- `import Enzyme`

### Frozen-ghost fix (applied to all 3 example files)

- `compute_f!` now calls `fill_halo_regions!(state)` before computing tendencies, ensuring the top
  halo `T[Nz+1] = 2·T_bc − T[Nz]` tracks the current interior temperature.

## Testing and verification

- Forward vs backward Euler: relative error ~1.7% after 10 timesteps with Δt = 300s (expected for
  implicit method with moderate Δt).
- Energy range evolves correctly: top boundary warms from prescribed surface temperature.

## Known limitations

- Dense Jacobian: O(N_z²) storage and O(N_z³) LU factorization per Newton iteration. Fine for CPU
  testing; the tridiagonal structure should be exploited for production (Phase B: linearized tridiagonal
  path via Oceananigans `BatchedTridiagonalSolver`).
- CTS cache `j` field must be a `DenseMatrix` for `lu(j)` — cannot store Jacobian as a Field directly.
  Fields are used for computation, flat Matrix for the CTS interface.
- No GPU support yet (CTS `lu` path is CPU-only).

## Future work

- Phase B: Linearized tridiagonal path via Oceananigans `BatchedTridiagonalSolver` (analytic Jacobian
  diagonals from the freeze curve and conductivity form).
- Differentiate `run_timesteps!` with Enzyme reverse-mode for full model sensitivities.
- GPU/Reactant compatibility.