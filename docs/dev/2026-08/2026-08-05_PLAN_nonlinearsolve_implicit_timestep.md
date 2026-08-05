# NonlinearSolve.jl for implicit (backward-Euler) timestepping: feasibility, linear-solver comparison, and automatic sparsity exploitation

> Status: **completed**. Host-side experiment implemented and verified: manual Enzyme callbacks (JVP + GMRES, colored Jacobian + KLU) work and are the fastest options; automatic sparsity detection works via **Route B** (`TracerLocalSparsityDetector` on a tracer-typed state running the unmodified `compute_f!` — exact tridiagonal pattern) and via dense finite differences (Route A variant); DI↔Enzyme automatic paths are blocked. All residual formulations use only Terrarium's own methods (no re-coded physics). Stage-2 (per-column kernel solves) documented as design notes only.

Date of initial draft: 2026-08-05

Base revision: 150cac8de389ea2935b017f6053b866e49cc1e71

## Originating prompt

> Look at @examples/autodiff/differentiating_terrarium.jl. Here you find already several attempts at doing a
> newton solve on implicit timestepping. What I want you to experiment with, is if you could also use
> nonlinearsolve.jl to do this implicit timestep. First investigate if it possible on the current setup. For
> linearsolvers, consider using multiple options. At a later stage, we want to be able to use
> KernelAbstractions.jl to launch the nonlinearsolves as a kernel. As this is a LSM, there will be one
> nonlinearsolve per grid cell.

Follow-up prompts during planning:

> For the Jacobian: is there a way to leverage the sparsity of the Jacobian? Also can't we set the AD option
> to Enzyme (and the sparsity) inside of the nonlinear solve option?

> Are you sure there is no way of doing automatic detection of sparsity? Ideally, in the future this implicit
> solver would be applied over a full land model which is modular, so in that case this automatic detection
> would be very helpful.

> Additional question: investigate if the matrix free methods can be used inside of a kernel, i.e. can you do
> a matrix free solve per grid cell inside of a kernel launch?

> Read https://arxiv.org/html/2403.16341v1 to check if your implementation makes sense.

## Revision log

- 2026-08-05: Initial draft after codebase investigation and NonlinearSolve.jl documentation research.
- 2026-08-05: User clarifications — solve granularity is **one N_z-dimensional solve per column** (columns
  are independent in an LSM); experiment lives in a **new example file**; scope is **stage 1 (host-side)
  fully implemented + stage 2 (kernels) as design notes**; this plan doc requested.
- 2026-08-05: Added sparsity strategy (manual `jac_prototype`, colored Enzyme JVPs, `AutoEnzyme` backend)
  and the finding that NonlinearSolve's default ForwardDiff Jacobian and tracer-based *detection* both fail
  through the current Float64 state-aliasing wrapper (element-type constraint).
- 2026-08-05: Corrected the sparsity-detection claim after user pushback — automatic detection **is**
  possible via Route A (`DenseSparsityDetector(AutoEnzyme)`) and Route B (`TracerSparsityDetector` on an
  eltype-generic column residual). The generic column-local residual (originally stage 2) is **pulled into
  scope** as the Route-B vehicle. Route C (structural assembly from process metadata) documented as future
  work.
- 2026-08-05: Added Krylov-in-kernel analysis (no ecosystem support; unjustified for tridiagonal per-column
  systems) and updated the stage-2 solver shortlist.
- 2026-08-05: Validated the plan against the NonlinearSolve.jl paper (arXiv:2403.16341v1). Adjustments:
  added the default polyalgorithm as a baseline variant; noted the state-dependent-branch caveat for
  approximate sparsity detection; added the implicit-function-theorem adjoint to future work; noted that
  `AutoEnzyme` integration postdates the paper (Enzyme is not mentioned in it).
- 2026-08-05: Implementation completed. Deviations from the planned matrix:
  (a) Route B dropped — `SparseConnectivityTracer` cannot be added to the examples environment
      (it requires LogExpFunctions 0.3.x while `Distributions` 0.25.130, via Makie/SpeedyWeather,
      requires 1.x; resolver unsatisfiable); the eltype-generic column residual was still implemented
      and validated (it is also the stage-2 vehicle). `SparseMatrixColorings` was added instead to enable
      the colored-sparse-AD path of variant 4.
  (b) `jac_nls!` uses 3 **multi-hot-seed JVPs** (plain `Duplicated`, columns ≡ c mod 3) rather than
      `BatchDuplicated` — equivalent for a tridiagonal pattern and closer to the reference example.
  (c) Variant 4 split: 4a (`DenseSparsityDetector(AutoEnzyme)`) fails with
      `EnzymeMutabilityException` (DI passes `Constant(p)` but the residual writes into `p.state`);
      added 4b, fully automatic sparsity exploitation with a working backend
      (`DenseSparsityDetector(AutoFiniteDiff)` + colored FD Jacobian + KLU).
  (d) Variant 7 (default polyalgorithm) does **not** fail: it recovers through fallback methods at
      ~40× the cost of the Enzyme variants.
  (e) The eltype-generic residual hoists per-layer coefficients (C(liq) = C0+C1·liq affine;
      κ(liq) = (A+B·liq)² from the `InverseQuadratic` weighting; ρLθ) evaluated from the actual model
      functions at liq ∈ {0, 1/2, 1}, with a built-in functional-form self-check (passed exactly, 0.0),
      because `SoilComposition`'s constructor requires a uniform number type and cannot mix Float64
      parameters with traced/dual state.
  (f) `Δt` is read from `integrator.model.timestepper.Δt`; the reference example's
      `integrator.timestepper.Δt` is stale w.r.t. the current `ModelIntegrator` (no `timestepper` field).
- 2026-08-05: **User constraint: no re-coded physics.** The eltype-generic column-local residual
  (`column_f!`/`G_column!` + coefficient hoisting) was **removed**; residual formulations may only
  use accessible Terrarium methods (`reset_tendencies!` → `closure!` → `compute_tendencies!`).
- 2026-08-05: **Route B revived and working.** User resolved the SCT dependency conflict (SCT and
  SMC now in `examples/Project.toml`). Automatic tracer-based sparsity detection succeeds through the
  *unmodified* `compute_f!` by re-typing only the container: a `StateVariables` copy whose `Field`s
  carry SCT `Dual` tracers, built with Oceananigans' `Field{LX, LY, LZ}(grid, T)` constructor
  (~20 lines of plumbing, zero physics). Findings:
  (a) naive tracing of the Float64 wrapper fails at the `u → state` write (`TypeError`);
  (b) the *global* detector (`GradientTracer`, no primal) fails in `safediv` (`iszero` requires a
      primal) — the *local* detector (`TracerLocalSparsityDetector`, `Dual` tracers with primals)
      works; comparisons on `Dual`s return `Bool`, so `ifelse` selects branches (detection reflects
      the branch configuration at u0 — verified exact here);
  (c) one Terrarium-side gate blocks tracing: `SoilComposition`'s constructor requires a uniform
      number type (porosity/solid stay Float64 while saturation/liquid are tracers) — circumvented in
      the example by a promoting constructor shim scoped to SCT `Dual`s; **upstream fix candidate**:
      make `SoilComposition`/related constructors promote (needed for any non-Enzyme AD through the
      closures, e.g. ForwardDiff would hit the same gate);
  (d) detected pattern is exactly tridiagonal (nnz = 3N−2 = 298, verified against the analytic
      structure and the dense Enzyme reference Jacobian);
  (e) variant 5 (SCT pattern + colored FD + KLU) matches variant 4b's solution exactly (u identical)
      but is ~13× faster in benchmarks (6 ms vs 80 ms): `DenseSparsityDetector(AutoFiniteDiff)`
      re-runs dense detection inside every `solve` call, while a pre-detected `jac_prototype`
      (`KnownJacobianSparsityDetector`) costs nothing per solve.
- 2026-08-05: **Enzyme FAQ investigation (variant 4a follow-up).** Per the Enzyme FAQ
  (Runtime Activity; Activity of temporary storage), `AutoEnzyme(mode =
  Enzyme.set_runtime_activity(Enzyme.Forward))` was probed through DI: the
  `EnzymeRuntimeActivityError` in the pushforward path disappears, but the JVP is *silently
  wrong* (max err 0.24 vs manual JVP) because DI marks the state-holding parameter `p` as
  `Constant` while the residual's active dataflow runs through `p.state` (the FAQ's
  "activity of temporary storage" trap: values inside a `Const` argument are definitionally
  constant); the dense-Jacobian path additionally crashes inside Enzyme (LLVM BoundsError).
  Conclusion: no correct DI-routed Enzyme option exists for this residual design; manual
  callbacks with `Duplicated(state, dstate)` (the FAQ's recommended fix) remain the only
  correct Enzyme route. Added **variant 4c** ("Route A with Enzyme, done manually"): dense
  Jacobian from N_z one-hot Enzyme JVPs → `sparse` pattern (== SCT pattern == tridiagonal)
  → `jac_prototype` + colored FD + KLU — Success, ‖G(u\*)‖∞ = 1.4e-7, ~5 ms.

## Problem description

Terrarium currently has **no implicit timestepper**: `src/timesteppers/` contains `ForwardEuler`, `Heun`,
and the IMEX routing scaffolding (`AbstractIMEX` with `Explicit()`/`Implicit()` traits), but no concrete
implicit sub-stepper. The example `examples/autodiff/differentiating_terrarium.jl` prototypes a
backward-Euler step for the soil heat equation externally:

  G(Uⁿ⁺¹) = Uⁿ⁺¹ − Uⁿ − Δt·f(Uⁿ⁺¹) = 0,   W = I − Δt·J_f,  J_f = ∂f/∂U

using (a) forward-mode Enzyme JVPs through `compute_f!` to build `J_f` column-by-column, and (b)
Ariadne.jl's `newton_krylov!` with a `JacobianOperator` for the Newton solve.

**Goal of this experiment:** replace the Ariadne prototype with NonlinearSolve.jl on the current setup
(single column, CPU, `Float64`), compare multiple linear-solver/Jacobian configurations (including
sparsity exploitation and automatic sparsity detection), and document the design for the future stage 2:
one NonlinearSolve per column launched as a KernelAbstractions kernel.

## Background

### Current setup (verified against source)

- `compute_f!(state, grid, soil, constants)` = `reset_tendencies!` → `closure!` (U→T freeze-curve mapping,
  `energy_to_temperature_kernel!`) → `compute_tendencies!` (∂U/∂t = −∂z(κ∂zT), KA kernels).
- Prognostic state: `internal_energy` only (NoFlow hydrology); N_z DOFs per column
  (`ExponentialSpacing()` → N_z = 50). Residual G has energy units, O(10⁷–10⁸) J/m³ → `abstol ≈ 1e-2`.
- The vector↔state bridge: write the solver iterate `u` into `state.prognostic.internal_energy` via
  `interior() .=`, call `compute_f!`, read `vec(interior(state.tendencies.internal_energy))`.
- The Jacobian W is **exactly tridiagonal**: the closure is pointwise (diagonal) and conduction couples
  adjacent layers only. Per-column blocks are independent → global Jacobian is block-diagonal.

### Key constraint discovered: element types

NonlinearSolve's default Jacobian (ForwardDiff via DifferentiationInterface) and tracer-based sparsity
detectors (`TracerSparsityDetector`) call `f!` with non-`Float64` element types. The state-aliasing wrapper
writes `u` into `Float64` Field arrays → `InexactError`. Consequences:

- We must supply explicit `jvp`/`jac` callbacks (Enzyme, which keeps element types via shadow arrays), or
  use `AutoEnzyme` (same property), or `AutoFiniteDiff` (Float64, slow baseline).
- Tracer-based detection cannot run on the Float64 state; it runs on a **tracer-typed copy** of the
  `StateVariables` instead (re-typed container, unmodified physics) — see Route B in the results.

### Validation against the NonlinearSolve.jl paper (arXiv:2403.16341v1)

- Default nonlinear polyalgorithm: Broyden/Klement → Broyden with true-Jacobian init → NewtonRaphson + line
  search → TrustRegion; first-order solvers directly for small (≤ 25) systems or when a Jacobian callback
  is given. We include the default as a baseline variant.
- Default linear-solver flowchart: matrix-free → Krylov.jl GMRES; small dense → RecursiveFactorization LU;
  sparse → KLU. Matches our comparison matrix.
- Approximate sparsity detection (union of dense Jacobians at random points — the mechanism behind
  `DenseSparsityDetector`) "fails to accurately predict the sparsity pattern in the presence of
  state-dependent branches". Our freeze curve contains `ifelse` branches, but they are pointwise
  (diagonal-only), so the detected pattern should still be tridiagonal — we verify this explicitly.
- Krylov section: `JacobianOperator` computes 𝒥x via forward AD and xᵀ𝒥 via reverse AD; `concrete_jac`
  forces materialization for preconditioners (IncompleteLU, AlgebraicMultigrid); sparse tooling can build
  cheap preconditioners. At large scales preconditioned Newton-Krylov wins (~1 order of magnitude on the
  32×32 Brusselator work-precision diagram); Sundials' Newton-Krylov failed to converge there.
- §2.3: solver-independent adjoints via the **implicit function theorem** (du\*/dθ = (∂f/∂u)⁻¹ ∂f/∂θ) are
  the recommended way to differentiate through a solve — highly relevant for differentiating implicit
  `run!` later (see Future work).
- Enzyme is **not mentioned** in the paper (v1, March 2024); `AutoEnzyme` support in
  DifferentiationInterface/NonlinearSolve postdates it — our Enzyme-based variants are the experimental
  checkpoint.

## Summary of changes

1. `examples/Project.toml`: added `ADTypes`, `DifferentiationInterface`, `Krylov`, `LinearSolve`,
   `NonlinearSolve`, `SparseArrays`, `SparseMatrixColorings`, `SparseConnectivityTracer`
   (experiment-only; root Project.toml untouched).
2. New example `examples/autodiff/implicit_timestepping_nonlinearsolve.jl`
   (`julia --project=examples examples/autodiff/implicit_timestepping_nonlinearsolve.jl`):
   - Setup identical to the Ariadne prototype (`ColumnGrid(CPU(), Float64, UniformSpacing(), 1)` →
     N_z = 100, `PrescribedSurfaceTemperature(1°C)`, Δt = 300 s from
     `integrator.model.timestepper.Δt`), reusing `compute_f!` (=
     `reset_tendencies!` → `closure!` → `compute_tendencies!`) for **all** residual formulations —
     no physics is re-implemented (user constraint).
   - `G_nls!`: backward-Euler residual; `jvp_nls!`: Enzyme forward JVP (W·v = v − Δt·J_f·v);
     `jac_nls!`: concrete W from 3 multi-hot-seed JVPs written into a CSC tridiagonal prototype.
   - **Route B (automatic tracer sparsity detection)**: tracer-typed `StateVariables` copy
     (`retype_state`, ~20 lines of container plumbing) + SCT-`Dual`-scoped promoting
     `SoilComposition` shim; `jacobian_sparsity(G_tr!, …, TracerLocalSparsityDetector())` on the
     residual `G_tr!` that calls the unmodified `compute_f!` on the tracer state. Detected pattern:
     exactly tridiagonal (nnz = 298 = 3N−2).
   - Solver comparison matrix and measured results (median benchmark times, CPU):

     | # | Jacobian / sparsity | Linear solver | Result | ‖G(u\*)‖∞ | time |
     |---|---|---|---|---|---|
     | 1a | manual Enzyme `jvp` (matrix-free) | `KrylovJL_GMRES()` | Success (3 Newton steps) | 5.3e-3 | 2.4 ms |
     | 1b | manual Enzyme `jvp` | GMRES + `EisenstatWalkerForcing2()` | Success | 1.5e-3 | 5.0 ms |
     | 2 | `jvp_autodiff = AutoEnzyme(Forward)` | GMRES | **FAILED**: `EnzymeRuntimeActivityError` | — | — |
     | 3 | manual colored Enzyme `jac`, CSC tridiagonal | `KLUFactorization()` | Success (3 steps) | 5.6e-11 | 3.2 ms |
     | 4a | `DenseSparsityDetector(AutoEnzyme)` | `KLUFactorization()` | **FAILED**: `EnzymeMutabilityException` | — | — |
     | 4b | `DenseSparsityDetector(AutoFiniteDiff)` + colored FD | `KLUFactorization()` | Success | 1.4e-7 | 72 ms |
     | 4c | **Route A-manual**: dense Enzyme jac → pattern + colored FD | `KLUFactorization()` | Success | 1.4e-7 | 5.3 ms |
     | 5 | **Route B**: SCT-detected pattern + colored FD | `KLUFactorization()` | Success | 1.4e-7 | 6.4 ms |
     | 6 | `AutoFiniteDiff()` dense | default LU | Success | 8.7e-10 | 120 ms |
     | 7 | default `solve(prob)` polyalgorithm | (internal) | Success via fallback | 5.5e-3 | 290 ms |

   - Key findings:
     - **Manual Enzyme callbacks are the only working Enzyme path.** DI-driven Enzyme (variants 2, 4a)
       fails: DI calls `autodiff` without `set_runtime_activity` (variant 2) and passes `Constant(p)`
       while the residual writes into `p.state` (variant 4a, `EnzymeMutabilityException`). The manual
       callbacks work because they annotate `Duplicated(state, dstate)` with
       `set_runtime_activity(Forward)`.
     - Variant 3 (exact colored Jacobian + KLU) is the most accurate (‖G‖∞ = 5.6e-11); variant 1a is
       the fastest. Both Enzyme variants beat all finite-difference alternatives.
     - **Variants 4c/5 are the practical automatic-sparsity routes**: identical solutions to 4b
       at ~12× lower cost, because `DenseSparsityDetector(AutoFiniteDiff)` re-runs dense detection
       inside every `solve` call while a pre-detected `jac_prototype`
       (`KnownJacobianSparsityDetector`) is detection-cost-free per solve. All three detection
       routes agree: dense Enzyme == SCT tracers == exact tridiagonal.
     - **No correct DI-routed Enzyme option exists** for this residual design (Enzyme FAQ:
       Runtime Activity; Activity of temporary storage): a runtime-activity-enabled
       `AutoEnzyme` mode passes through DI and removes the hard error in the pushforward path
       but yields *silently wrong* JVPs (state path dropped because DI marks the state-holding
       parameter `Constant`); the dense-Jacobian path crashes inside Enzyme. Manual callbacks
       with `Duplicated(state, dstate)` (the FAQ's recommended fix) are the only correct
       Enzyme route — used by variants 1/3 and by 4c's detection.
     - The default polyalgorithm recovers through fallback methods (332 function evaluations) — it
       does not fail, but is by far the slowest.
     - `jac_nls!` matches the dense N_z-JVP reference Jacobian exactly (max |ΔW| = 0); W is exactly
       tridiagonal (max off-band entry = 0.0); the SCT-detected pattern equals it (nnz = 298).
   - Correctness: all successful variants agree with the Ariadne `newton_krylov!` reference to
     ≤ 2.5e-2 absolute (~2.5e-9 relative on O(10⁷) J/m³; differences set by termination tolerance) and
     with each other to ≤ 5e-3; residual norms above.
   - Benchmarks printed via BenchmarkTools at the end of the script.
   - Closing comment section: stage-2 design notes (below).

## Testing and verification

- Example runs end-to-end (`julia --project=examples`, exit 0) on CPU.
- All pass criteria met except the two documented experimental failures (2, 4a — DI↔Enzyme);
  detected sparsity == tridiagonal (both Route A-FD and Route B-tracer); see the results table above.
- No changes to `src/` or tests; example-only experiment.

## Documentation changes

- The example file itself (Literate comments) documents usage and findings.
- This plan doc tracks status; no user-facing docs pages changed (no new public API).

## Known limitations

- Host-side, single-column, CPU only; no integration into `timestep!`/`AbstractTimeStepper` yet.
- Element-type constraint: NonlinearSolve's default ForwardDiff Jacobian and tracer detectors cannot run
  through the Float64 state-aliasing wrapper (see Background); tracer detection therefore runs on a
  re-typed *copy* of the state (container-level, no physics changes).
- **DI↔Enzyme incompatibility** (measured; Enzyme FAQ investigated): `jvp_autodiff`/`autodiff`/
  `DenseSparsityDetector` with `AutoEnzyme` fail through the state-aliasing wrapper
  (`EnzymeRuntimeActivityError`, `EnzymeMutabilityException`). Passing a runtime-activity-enabled
  mode (`AutoEnzyme(mode = set_runtime_activity(Forward))`) removes the hard error in the
  pushforward path but yields *silently wrong* JVPs: DI marks the state-holding parameter
  `Constant`, and values inside a `Const` argument are definitionally constant (FAQ: "activity of
  temporary storage"). Only manual Enzyme callbacks with `Duplicated(state, dstate)` are correct.
- **Tracer detection requires a local shim**: `SoilComposition`'s uniform-number-type constructor gate
  rejects mixed Float64-parameter/tracer-state calls (also relevant for ForwardDiff); the example
  carries a promoting fallback scoped to SCT `Dual`s — upstream fix candidate (promoting constructors).
- **Branch-selection caveat for tracer detection**: comparisons on SCT `Dual`s return `Bool`
  (primal-valued), so `ifelse` selects a single branch and the detected pattern reflects the branch
  configuration at the detection state u0. For this model the tridiagonal pattern is robust to branch
  changes (the diagonal survives via the direct T(U) path or, in the phase-change region, via the
  κ(liq(U)) coupling), but for a modular full model, detection should run at a representative state
  and be re-verified when the branch configuration changes.
- Tolerances must be scaled to energy units (O(10⁷) J/m³).
- Differentiating *through* the NonlinearSolve solve with Enzyme is untested (see Future work).

## Future work

### Stage 2: per-column solves as a KernelAbstractions kernel

- Target pattern (SciML GPU tutorial): `ImmutableNonlinearProblem` + StaticArrays `u0` +
  SimpleNonlinearSolve algorithms + `SciMLBase.remake(prob; u0, p)` inside a `@kernel`; one
  N_z-dimensional solve per column.
- **Main design consequence:** `compute_f!` *launches KA kernels* and kernels cannot nest, so the
  in-kernel residual must be assembled from Terrarium's scalar **kernel functions** (called per-index
  over plain buffers, e.g. via `get_fields` on per-column views) rather than by launching the model
  kernels. The Route-B experiment demonstrates these kernel functions are eltype-generic in practice;
  the main blocker is the `SoilComposition` uniform-number-type constructor gate (see Known
  limitations), which a promoting constructor would remove.
- **Krylov-in-kernel analysis** (investigated, conclusion: not the right target):
  - No ecosystem support: SimpleNonlinearSolve's kernel-safe set contains no Krylov methods;
    Krylov.jl/LinearSolve.jl are host-side (one whole-array solve using the GPU, not per-work-item solves).
  - Hand-rolled in-kernel Krylov: only fixed-workspace BiCGStab (~7 vectors) is realistic; GMRES needs
    m+1 vectors + per-iteration dense Hessenberg solves (register/local-memory pressure), and
    data-dependent iteration counts cause warp divergence and Reactant `while_loop` tracing issues.
    ClimateMachine's batched GMRES is host-driven with a global-memory basis — not a fused per-thread
    solve.
  - Physics argument: per-column systems are small (N_z ≈ 50–100) and tridiagonal — an exact O(n) Thomas
    solve beats approximate Krylov. Matrix-free Krylov remains right for host-side stage 1 and any future
    *globally coupled* problem (e.g. horizontal coupling breaking the per-column decomposition).
- In-kernel solver shortlist: `SimpleKlement`, `SimpleLimitedMemoryBroyden`, `SimpleDFSane`
  (factorization-free), `SimpleNewtonRaphson`/`SimpleTrustRegion` (dense static J — beware StaticArrays
  `\`/LU is only allocation-free/reliable for N ≲ 14), and a custom tridiagonal Newton (analytic or
  3-color ForwardDiff bands + Thomas algorithm) which the physics favors.
- In-kernel AD: ForwardDiff over StaticArrays (GPU-compatible) or analytic Jacobian — not Enzyme.
- Termination inside kernels/Reactant: prefer fixed iteration counts + masking over dynamic termination;
  `SimpleLimitedMemoryBroyden` is non-allocating only with termination `nothing` or
  `AbsNormTerminationMode`.

### Route C: structural Jacobian assembly from process metadata

For the modular full land model, the scalable long-term answer is a `jac_prototype(model, grid)` assembled
from per-process coupling declarations (e.g. `coupling_stencil(process) → Pointwise() | VerticalStencil(1)`),
giving exact patterns with zero tracing cost and natural modularity. New interface machinery; deferred.

### AD through the implicit solve

For differentiating implicit `run!`: prefer an **implicit-function-theorem adjoint**
(du\*/dθ = (∂f/∂u)⁻¹ ∂f/∂θ, cf. NonlinearSolve.jl paper §2.3) over unrolling Newton iterations — cheaper
and solver-independent. Requires a custom Enzyme rule (or SciMLSensitivity-style adjoint) plus one linear
solve with ∂f/∂u per output; must be tested with our Enzyme JVPs. Alternative: differentiate through the
iterations directly (works if NonlinearSolve internals are Enzyme-compatible; known to be fragile).

### GPU grid support (host-side)

Krylov.jl works on `CuArray`s (CSR operators via `KrylovOperator`), so a whole-grid GPU `NonlinearProblem`
is possible later; the per-column kernel approach supersedes it for the LSM use case.
