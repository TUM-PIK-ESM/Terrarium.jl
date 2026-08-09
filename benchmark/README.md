# Benchmarks

Performance benchmarks for Terrarium.jl, collected across architectures. Each architecture's results live in its own section below; the overview table compares the resolution sweeps across all architectures benchmarked so far.

Every configuration is run for a fixed number of time steps without output. Initialization and — under Reactant — XLA compilation happen before the clock starts and are reported separately; the timed section covers only the stepping loop. A run takes a few seconds and is timed once, so it measures the mean rather than the minimum: deviations of ±20% between repetitions are normal, and `quick` mode — which shortens every run — is noisier still, so treat its numbers as a smoke test rather than a measurement. Use the `long` mode for numbers worth quoting.

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
| land | 128 | 158 |
| land | 512 | 38 |
| land | 2048 | 11 |
| land | 4608 | 2.2 |
| land | 8192 | 2.87 |
| land | 18432 | — |
| land | 41472 | — |
| land | 73728 | — |
| land | 165888 | — |
| land_no_vegetation | 128 | 181 |
| land_no_vegetation | 512 | 37 |
| land_no_vegetation | 2048 | 8.72 |
| land_no_vegetation | 4608 | 3.91 |
| land_no_vegetation | 8192 | 2.6 |
| land_no_vegetation | 18432 | — |
| land_no_vegetation | 41472 | — |
| land_no_vegetation | 73728 | — |
| land_no_vegetation | 165888 | — |
| soil_heat | 128 | 269 |
| soil_heat | 512 | 97 |
| soil_heat | 2048 | 26 |
| soil_heat | 4608 | 12 |
| soil_heat | 8192 | 6.97 |
| soil_heat | 18432 | — |
| soil_heat | 41472 | — |
| soil_heat | 73728 | — |
| soil_heat | 165888 | — |

## Architecture: `cpu-arm`

Created for Terrarium.jl v0.1.4 on Mon, 03 Aug 2026 17:44:34 in `quick` mode (0.25x time steps, 1 thread(s)).

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
| land | 3.75° | 4608 | 600 | 54 | 2.52 | 651 | 0.212 | 6.85 MiB |
| land_no_vegetation | 3.75° | 4608 | 600 | 54 | 5.85 | 281 | 0.492 | 5.76 MiB |
| soil_heat | 3.75° | 4608 | 600 | 54 | 14 | 119 | 1.17 | 3.85 MiB |

### Land model, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 22.50° | 128 | 600 | 500 | 158 | 10.4 | 0.369 | 203.62 KiB | ok |
| land | 11.25° | 512 | 600 | 488 | 38 | 43.7 | 0.352 | 787.12 KiB | ok |
| land | 5.62° | 2048 | 600 | 122 | 11 | 152 | 0.405 | 3.05 MiB | ok |
| land | 3.75° | 4608 | 600 | 54 | 2.2 | 748 | 0.185 | 6.85 MiB | ok |
| land | 2.81° | 8192 | 600 | 30 | 2.87 | 573 | 0.429 | 12.17 MiB | ok |
| land | 1.88° | 18432 | — | — | — | — | — | — | skipped |
| land | 1.25° | 41472 | — | — | — | — | — | — | skipped |
| land | 0.94° | 73728 | — | — | — | — | — | — | skipped |
| land | 0.62° | 165888 | — | — | — | — | — | — | skipped |

### Land model without vegetation, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land_no_vegetation | 22.50° | 128 | 600 | 500 | 181 | 9.05 | 0.424 | 171.16 KiB | ok |
| land_no_vegetation | 11.25° | 512 | 600 | 488 | 37 | 44.9 | 0.342 | 661.66 KiB | ok |
| land_no_vegetation | 5.62° | 2048 | 600 | 122 | 8.72 | 188 | 0.326 | 2.56 MiB | ok |
| land_no_vegetation | 3.75° | 4608 | 600 | 54 | 3.91 | 420 | 0.329 | 5.76 MiB | ok |
| land_no_vegetation | 2.81° | 8192 | 600 | 30 | 2.6 | 632 | 0.389 | 10.23 MiB | ok |
| land_no_vegetation | 1.88° | 18432 | — | — | — | — | — | — | skipped |
| land_no_vegetation | 1.25° | 41472 | — | — | — | — | — | — | skipped |
| land_no_vegetation | 0.94° | 73728 | — | — | — | — | — | — | skipped |
| land_no_vegetation | 0.62° | 165888 | — | — | — | — | — | — | skipped |

### Soil heat conduction, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| soil_heat | 22.50° | 128 | 600 | 500 | 269 | 6.11 | 0.629 | 114.63 KiB | ok |
| soil_heat | 11.25° | 512 | 600 | 488 | 97 | 16.9 | 0.908 | 443.13 KiB | ok |
| soil_heat | 5.62° | 2048 | 600 | 122 | 26 | 62.4 | 0.984 | 1.72 MiB | ok |
| soil_heat | 3.75° | 4608 | 600 | 54 | 12 | 135 | 1.02 | 3.85 MiB | ok |
| soil_heat | 2.81° | 8192 | 600 | 30 | 6.97 | 236 | 1.04 | 6.85 MiB | ok |
| soil_heat | 1.88° | 18432 | — | — | — | — | — | — | skipped |
| soil_heat | 1.25° | 41472 | — | — | — | — | — | — | skipped |
| soil_heat | 0.94° | 73728 | — | — | — | — | — | — | skipped |
| soil_heat | 0.62° | 165888 | — | — | — | — | — | — | skipped |

### Number of soil layers

| Configuration | Res | Columns | L | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | 10 | 600 | 163 | 12 | 132 | 0.348 | 3.68 MiB |
| land | 3.75° | 4608 | 20 | 600 | 82 | 7.52 | 218 | 0.422 | 5.26 MiB |
| land | 3.75° | 4608 | 30 | 600 | 54 | 4.94 | 333 | 0.415 | 6.85 MiB |
| land | 3.75° | 4608 | 60 | 600 | 27 | 2.54 | 646 | 0.428 | 11.60 MiB |
| land | 3.75° | 4608 | 100 | 600 | 16 | 1.51 | 1090 | 0.423 | 17.94 MiB |

### Number format, Float32 vs Float64

| Configuration | NF | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | Float32 | 3.75° | 4608 | 600 | 54 | 4.38 | 375 | 0.368 | 6.85 MiB |
| land | Float64 | 3.75° | 4608 | 600 | 54 | 3.49 | 471 | 0.293 | 13.69 MiB |

### Time stepper

| Configuration | Variant | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | ForwardEuler | 3.75° | 4608 | 600 | 54 | 5.18 | 317 | 0.436 | 6.85 MiB |
| land | Heun | 3.75° | 4608 | 600 | 54 | 1.09 | 1510 | 0.0915 | 6.85 MiB |
