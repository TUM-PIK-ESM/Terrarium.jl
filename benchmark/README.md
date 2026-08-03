# Benchmarks

Performance benchmarks for Terrarium.jl, collected across architectures. Each architecture's results live in its own section below; the overview table compares the resolution sweeps across all architectures benchmarked so far.

Every configuration is run for a fixed number of time steps without output. Initialization and — under Reactant — XLA compilation happen before the clock starts and are reported separately; the timed section covers only the stepping loop. A run takes a few seconds and is timed once, so it measures the mean rather than the minimum: deviations of ±20% between repetitions are normal. Use the `long` mode for numbers that need to be more robust than that.

### Explanation

- Configuration: model setup, see `model_configurations.jl`
- NF: number format, default `Float32`
- Res: approximate horizontal resolution in degrees latitude
- Columns: number of land columns; the benchmark grids are full Gaussian grids with an all-land mask
- L: number of soil layers, default 30
- Δt: time step (s)
- SYPD: simulated years per wallclock day
- ms/step: wallclock milliseconds per time step
- Mcell-steps/s: millions of (column × layer) updates per second — the resolution-independent throughput
- Memory: state and tendency fields, including halos
- Compile: one-off XLA compilation time of the stepping loop (Reactant only)
- Status: `ok`, `skipped` (above the mode's resolution cap), `unstable` (state went non-finite) or `failed`

### Running the benchmarks

From `benchmark/`:

```
julia --project=. manual_benchmarking.jl                # CPU (auto-labelled cpu-arm or cpu-x86)
julia --project=. manual_benchmarking.jl gpu            # CUDA GPU
julia --project=. manual_benchmarking.jl reactant-cpu   # Reactant/XLA, CPU backend
julia --project=. manual_benchmarking.jl reactant-gpu   # Reactant/XLA, CUDA backend
```

A second argument controls the duration: `quick` (0.25x steps, sweeps capped at 8192 columns), `long` (10x steps), or a numeric timestep multiplier. Each run updates only its own architecture's section here; the other architectures are preserved in `assets/benchmark_results.json`.

## Overview: resolution sweeps across architectures

Simulated years per wallclock day (SYPD) for each model configuration across horizontal resolutions, one column per architecture. `—` means the architecture has not been benchmarked yet, skipped that resolution, or failed to run that configuration (see the per-architecture sections for which). Comparison figures are on the documentation's Benchmarks page.

| Configuration | Columns | cpu-arm |
| --- | --- | --- |
| land | 128 | 136 |
| land | 512 | 39 |
| land | 2048 | 11 |
| land | 4608 | 2.55 |
| land | 8192 | 2.6 |
| land | 18432 | — |
| land | 41472 | — |
| land | 73728 | — |
| land | 165888 | — |
| land_no_vegetation | 128 | 167 |
| land_no_vegetation | 512 | 46 |
| land_no_vegetation | 2048 | 11 |
| land_no_vegetation | 4608 | 5.03 |
| land_no_vegetation | 8192 | 2.82 |
| land_no_vegetation | 18432 | — |
| land_no_vegetation | 41472 | — |
| land_no_vegetation | 73728 | — |
| land_no_vegetation | 165888 | — |
| soil_heat | 128 | 403 |
| soil_heat | 512 | 104 |
| soil_heat | 2048 | 26 |
| soil_heat | 4608 | 11 |
| soil_heat | 8192 | 6.47 |
| soil_heat | 18432 | — |
| soil_heat | 41472 | — |
| soil_heat | 73728 | — |
| soil_heat | 165888 | — |

## Architecture: `cpu-arm`

Created for Terrarium.jl v0.1.4 on Mon, 03 Aug 2026 16:54:09 in `quick` mode (0.25x time steps, 1 thread(s)).

### Machine details

```julia
julia> versioninfo()
Julia Version 1.11.9
Commit 53a02c0720c (2026-02-06 00:27 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: macOS (arm64-apple-darwin24.0.0)
  CPU: 8 × Apple M3
  WORD_SIZE: 64
  LLVM: libLLVM-16.0.6 (ORCJIT, apple-m2)
Threads: 1 default, 0 interactive, 1 GC (on 4 virtual cores)
```


### Model configurations, default resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | 600 | 54 | 2.78 | 591.0 | 0.234 | 6.85 MiB |
| land_no_vegetation | 3.75° | 4608 | 600 | 54 | 5.45 | 301.0 | 0.459 | 5.76 MiB |
| soil_heat | 3.75° | 4608 | 600 | 54 | 13 | 128.0 | 1.08 | 3.85 MiB |

### Land model, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 22.50° | 128 | 600 | 500 | 136 | 12.1 | 0.318 | 203.62 KiB | ok |
| land | 11.25° | 512 | 600 | 488 | 39 | 41.6 | 0.369 | 787.12 KiB | ok |
| land | 5.62° | 2048 | 600 | 122 | 11 | 153.0 | 0.401 | 3.05 MiB | ok |
| land | 3.75° | 4608 | 600 | 54 | 2.55 | 645.0 | 0.214 | 6.85 MiB | ok |
| land | 2.81° | 8192 | 600 | 30 | 2.6 | 632.0 | 0.389 | 12.17 MiB | ok |
| land | 1.88° | 18432 | — | — | — | — | — | — | skipped |
| land | 1.25° | 41472 | — | — | — | — | — | — | skipped |
| land | 0.94° | 73728 | — | — | — | — | — | — | skipped |
| land | 0.62° | 165888 | — | — | — | — | — | — | skipped |

### Land model without vegetation, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land_no_vegetation | 22.50° | 128 | 600 | 500 | 167 | 9.85 | 0.39 | 171.16 KiB | ok |
| land_no_vegetation | 11.25° | 512 | 600 | 488 | 46 | 36.1 | 0.425 | 661.66 KiB | ok |
| land_no_vegetation | 5.62° | 2048 | 600 | 122 | 11 | 145.0 | 0.423 | 2.56 MiB | ok |
| land_no_vegetation | 3.75° | 4608 | 600 | 54 | 5.03 | 326.0 | 0.424 | 5.76 MiB | ok |
| land_no_vegetation | 2.81° | 8192 | 600 | 30 | 2.82 | 583.0 | 0.422 | 10.23 MiB | ok |
| land_no_vegetation | 1.88° | 18432 | — | — | — | — | — | — | skipped |
| land_no_vegetation | 1.25° | 41472 | — | — | — | — | — | — | skipped |
| land_no_vegetation | 0.94° | 73728 | — | — | — | — | — | — | skipped |
| land_no_vegetation | 0.62° | 165888 | — | — | — | — | — | — | skipped |

### Soil heat conduction, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| soil_heat | 22.50° | 128 | 600 | 500 | 403 | 4.07 | 0.943 | 114.63 KiB | ok |
| soil_heat | 11.25° | 512 | 600 | 488 | 104 | 15.8 | 0.971 | 443.13 KiB | ok |
| soil_heat | 5.62° | 2048 | 600 | 122 | 26 | 62.7 | 0.98 | 1.72 MiB | ok |
| soil_heat | 3.75° | 4608 | 600 | 54 | 11 | 148.0 | 0.937 | 3.85 MiB | ok |
| soil_heat | 2.81° | 8192 | 600 | 30 | 6.47 | 254.0 | 0.968 | 6.85 MiB | ok |
| soil_heat | 1.88° | 18432 | — | — | — | — | — | — | skipped |
| soil_heat | 1.25° | 41472 | — | — | — | — | — | — | skipped |
| soil_heat | 0.94° | 73728 | — | — | — | — | — | — | skipped |
| soil_heat | 0.62° | 165888 | — | — | — | — | — | — | skipped |

### Number of soil layers

| Configuration | Res | Columns | L | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | 10 | 600 | 163 | 12 | 134.0 | 0.344 | 3.68 MiB |
| land | 3.75° | 4608 | 20 | 600 | 82 | 6.6 | 249.0 | 0.37 | 5.26 MiB |
| land | 3.75° | 4608 | 30 | 600 | 54 | 4.59 | 358.0 | 0.386 | 6.85 MiB |
| land | 3.75° | 4608 | 60 | 600 | 27 | 2.42 | 678.0 | 0.408 | 11.60 MiB |
| land | 3.75° | 4608 | 100 | 600 | 16 | 1.44 | 1140.0 | 0.403 | 17.94 MiB |

### Number format, Float32 vs Float64

| Configuration | NF | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | Float32 | 3.75° | 4608 | 600 | 54 | 4.54 | 362.0 | 0.382 | 6.85 MiB |
| land | Float64 | 3.75° | 4608 | 600 | 54 | 3.17 | 519.0 | 0.267 | 13.69 MiB |

### Time stepper

| Configuration | Variant | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | ForwardEuler | 3.75° | 4608 | 600 | 54 | 4.69 | 350.0 | 0.395 | 6.85 MiB |
| land | Heun | 3.75° | 4608 | 600 | 54 | 1.16 | 1420.0 | 0.0974 | 6.85 MiB |

