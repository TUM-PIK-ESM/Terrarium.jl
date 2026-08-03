# Benchmark suite for Terrarium

> Status: **in progress**. Harness, model configuration registry, results store, README generator and
> documentation page are implemented; `cpu-arm` and `reactant-cpu` results collected locally, the
> remaining four architectures still need to be run on their respective machines.

Date of initial draft: 2026-08-03

Base revision: e015a29db9d97090edfd112a5eea165640134404

## Originating prompt

> We want to set up a benchmark suite for Terrarium.
> Inspect SpeedyWeather.jl's benchmark for an example, you can find them in the repository in the
> benchmark folder.
>
> * We need to run on ARM (normal), ARM (Reactant), x86 (normal), x86 (Reactant), GPU (normal),
>   GPU (Reactant)
> * Reactant needs to use `Reactant.set_default_backend("gpu")`
> * We need a list of model configurations to run. They should take in the number of gridpoints /
>   resolution
> * We want a similar documentation upload of the benchmarks

> Also introduce a modifier to run the tests longer for a more robust test, or short for a quick
> benchmark

## Revision log

- 2026-08-03: initial draft.
- 2026-08-03: model configuration list narrowed (in review) to the coupled `LandModel` with Richards
  hydrology and snow, plus a no-vegetation variant; `:soil_heat` retained only as a Reactant
  reference configuration. Resolution sweep extended to `nlat_half = 144` (~166k columns) on a
  synthetic all-land mask. The stale `test/benchmarks/gpu/` scripts are left untouched.

## Problem description

Terrarium has no reproducible performance benchmark. The only existing artefact,
`test/benchmarks/gpu/soil_heat_hydrology_global.jl`, is a one-off scaling script written against the
pre-`ModelIntegrator` API (`StateVariables(model)`, `timestep!(state, Δt)`) and no longer runs.

At the same time Terrarium has four execution paths that must each stay fast — plain CPU, CUDA GPU,
Reactant/XLA on CPU, Reactant/XLA on GPU — across two CPU ISAs (ARM and x86), i.e. six benchmark
targets. There is currently no way to state what any of them costs, no baseline against which a
performance regression could be noticed, and no published numbers for users deciding which backend
to run on.

## Background

SpeedyWeather.jl solves the same problem in `SpeedyWeather/benchmark`:

- `manual_benchmarking.jl <arch>` runs the suites for one architecture,
- merges the results into a shared `assets/benchmark_results.json` keyed by architecture label,
- regenerates `README.md` from the *whole* store, so a run on one machine never clobbers another's,
- and `docs/generate_benchmarks_page.jl` renders the same JSON into a `Benchmarks` documentation
  page with cross-architecture comparison figures at doc-build time.

That design is adopted here essentially unchanged, because the awkward part of this task — six
architectures that can never be measured on one machine at one time — is exactly what the JSON store
plus per-architecture regeneration solves.

## Summary of changes

### `benchmark/` (new)

| File | Purpose |
| --- | --- |
| `Project.toml` | Benchmark environment; `Terrarium` via `[sources]` path `..` |
| `model_configurations.jl` | `Val`-dispatched registry of resolution-parametrized model configurations |
| `benchmark_suite.jl` | `BenchmarkSuite` type, timing protocol, markdown writers |
| `define_benchmarks.jl` | The suites: which configurations at which resolutions |
| `manual_benchmarking.jl` | Driver: architecture/mode parsing, backend loading, JSON merge, README regeneration |
| `assets/benchmark_results.json` | Committed results store, one record per architecture |
| `README.md` | Generated: overview table + one section per architecture |

### Architecture matrix

Selected by the first CLI argument. The label is the JSON key, the README section title and the
documentation table column.

| Argument | Label | Architecture | Backends loaded |
| --- | --- | --- | --- |
| *(none)* / `cpu` | `cpu-arm` / `cpu-x86` (from `Sys.ARCH`) | `CPU()` | — |
| `gpu` | `gpu-nvidia` | `GPU()` | `CUDA` |
| `reactant-cpu` | `reactant-cpu` | `ReactantState()` | `CUDA`, `Reactant` + `set_default_backend("cpu")` |
| `reactant-gpu` | `reactant-gpu` | `ReactantState()` | `CUDA`, `Reactant` + `set_default_backend("gpu")` |

`CUDA` is loaded for *both* Reactant variants: it provides the KernelAbstractions↔Reactant glue and
is required even on the CPU backend. All backend packages are loaded before `using Terrarium` so the
package extensions register.

### Duration modifier

Second CLI argument, so one suite serves both a smoke check and a publication run:

```
julia --project=. manual_benchmarking.jl gpu quick   # 0.25x steps, sweep capped at 8192 columns
julia --project=. manual_benchmarking.jl gpu         # default, each timed run a few seconds
julia --project=. manual_benchmarking.jl gpu long    # 10x steps, full sweep
julia --project=. manual_benchmarking.jl gpu 3.5     # explicit timestep multiplier, full sweep
```

`BenchmarkMode` carries `timestep_multiplier` and `max_npoints`. The mode is recorded in the JSON
`meta` and printed in the README and documentation heading for each architecture: SYPD is a rate and
so is comparable across modes, but the noise level is not, and the two must not be silently mixed.

### Model configurations

`build_model(::Val{:name}, arch, NF; nlat_half, nz)` returns
`(; model, boundary_conditions, initializers, Δt)`, following the registry pattern already used by
`test/reactant/setup.jl`. `arch`, `nlat_half` and `nz` are the only things that vary between runs.

Horizontal resolution is `RingGrids.FullGaussianGrid(nlat_half)` with the default all-land mask
(`npoints = 8·nlat_half²`), so no asset download is needed on any benchmark machine — the same
choice `test/reactant/setup.jl` makes for CI. Vertical resolution is `ExponentialSpacing(N = nz)`.

- `:land` — fully coupled `LandModel`: `VegetationCarbon`, `SoilEnergyWaterCarbon` with
  `SoilHydrology(NF, RichardsEq())`, `SingleLayerSnow`, and the default surface energy balance,
  surface hydrology and prescribed atmosphere. `Δt = 600 s`; the coupled model raises a `DomainError`
  in the turbulent-flux thermodynamics for Δt ≳ 1800 s.
- `:land_no_vegetation` — identical but `vegetation = nothing`, which also selects
  `BareGroundEvaporation` and `NoCanopyInterception`. The pair isolates the cost of the
  vegetation/canopy path.
- `:soil_heat` — the minimal heat-conduction `SoilModel` from `test/reactant/setup.jl`. This is a
  *reference* configuration, not part of the headline comparison: it is the only configuration
  currently proven to compile under Reactant, so an empty Reactant column can be attributed to the
  land physics rather than to the harness.

### Suites

| Key | Title | Sweep |
| --- | --- | --- |
| `:bench100` | Model configurations, default resolution | all configurations, `nlat_half = 24`, `nz = 30` |
| `:bench200` | Land model, horizontal resolution | `:land`, `nlat_half ∈ [4, 8, 16, 24, 32, 48, 72, 96, 144]` — the cross-architecture overview |
| `:bench201` | Land model without vegetation, horizontal resolution | `:land_no_vegetation`, same sweep |
| `:bench300` | Vertical resolution | `:land`, `nz ∈ [10, 20, 30, 60, 100]` |
| `:bench400` | Number format | `:land`, `Float32` vs `Float64` |
| `:bench500` | Time stepper | `:land`, `ForwardEuler` / `Heun` / `IMEX` |

Configurations above the mode's `max_npoints` are skipped, logged, and rendered as `—`. Nothing is
dropped silently.

### Timing protocol

Per configuration: build and `initialize` (timed separately, reported as *init*, excluded from
SYPD); compute state memory by summing `length(parent(f)) · sizeof(eltype(f))` over the model's
state and tendency fields (`Base.summarysize` is not used — it undercounts device arrays by orders
of magnitude); warm-up run; re-`initialize!`; timed run with `time()` around it and an architecture
synchronization before the clock is stopped.

- CPU: no synchronization needed.
- GPU: `CUDA.synchronize()`.
- Reactant: `run!` recompiles on every call (its `compiled_run!` argument is never handed back to the
  caller), so the harness compiles once itself, exactly as `TerrariumReactantExt` does —
  `@compile raise = true raise_first = true sync = true run_timesteps!(integrator, Δt, nsteps, false)` —
  reports the compile time as its own metric, and then calls the compiled function twice (warm-up,
  timed). `sync = true` makes the timed call synchronous.

Step counts follow `nsteps = multiplier · clamp(round(Int, C / (npoints·nz)), 20, 2000)`, with `C`
chosen so a mid-size CPU configuration takes a few seconds.

Any failure of a single configuration is caught, logged with the exception, and recorded as `—`, so
one non-compiling configuration cannot cost an entire architecture its results.

### Metrics

SYPD (simulated years per wallclock day, `Δt·nsteps / (elapsed·365.25)` — the same definition
SpeedyWeather uses, so the two models' tables are directly comparable), ms/step, throughput in
cell-steps per second, state memory, initialization time, and (Reactant only) compile time.

### Documentation

`docs/generate_benchmarks_page.jl` is included from `docs/make.jl` before `makedocs` and writes
`docs/src/benchmarks.md` from the JSON store, together with CairoMakie figures (SYPD vs number of
columns, log-log, one line per architecture, one figure per configuration, with reference lines at
5°/2°/1°/0.5°). Section headings are suffixed with the architecture label so Documenter's `@ref`
resolver cannot confuse them with other pages. A missing or empty JSON store yields a placeholder
page, so the documentation still builds on a fresh checkout. `docs/Project.toml` gains `JSON3`;
`docs/make.jl` gains a top-level `"Benchmarks"` page.

## Testing and verification

Benchmarks are not part of CI (they are wallclock measurements on dedicated hardware). Verification
is manual:

1. `julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'`
2. `cd benchmark && julia --project=. manual_benchmarking.jl cpu quick` — populates the `cpu-arm`
   (or `cpu-x86`) section of `benchmark/README.md` and writes `assets/benchmark_results.json`.
3. `julia --project=. manual_benchmarking.jl reactant-cpu quick` — exercises the compile path.
4. `julia --project=docs docs/make.jl --draft` — renders `docs/src/benchmarks.md` and its figures.
5. Running twice must change only the record for the architecture just measured.

## Documentation changes

- New generated page `docs/src/benchmarks.md`, registered as a top-level `Benchmarks` entry.
- New `benchmark/README.md` (generated) describing the metric definitions and how to reproduce.
- `AGENTS.md` gains a short `## Benchmarks` section: how to run the suite, and the convention that
  the committed README numbers are the regression baseline (±20% is noise; larger deviations should
  be reported).

## Known limitations

- **Reactant coverage of the land model is unproven.** `:land` is expected to fail to compile on the
  static root-fraction path (`exp(::TracedRArray)` from the `FunctionField` `sum(dims = 3)` in
  `StaticExponentialRootDistribution`); `:land_no_vegetation` is untested. The harness records such
  failures as `—` rather than hiding them, and `:soil_heat` distinguishes a physics gap from a
  harness bug.
- **Single-threaded CPU numbers.** Thread scaling is not swept; runs inherit whatever
  `JULIA_NUM_THREADS` the machine provides, and the value is recorded in the machine info block.
- **Synthetic all-land mask.** Column counts are ~3x a real land-sea mask at the same `nlat_half`.
  This is deliberate (reproducibility, no asset download) but means absolute SYPD is pessimistic
  relative to a realistic land simulation.
- **Timings are means, not minima.** As in SpeedyWeather, a run of several seconds is timed once;
  ±20% variation between repetitions is normal. The `long` mode exists for when that is not good
  enough.

## Future work

- `test/benchmarks/gpu/soil_heat_hydrology_global.jl` and its `.sbatch` companion are written
  against the pre-`ModelIntegrator` API and no longer run. They are left untouched by this change;
  they should either be ported to the new harness or removed in a separate PR.
- Automated regression checking (compare against the committed JSON in CI on a fixed runner).
- A thread-scaling suite for the CPU architectures.
- Per-kernel benchmarks (analogous to SpeedyWeather's individual dynamics functions suite) once the
  full-model numbers are stable.
