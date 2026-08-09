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

| Configuration | Columns | cpu-x86 | gpu-nvidia | reactant-cpu | reactant-gpu |
| --- | --- | --- | --- | --- | --- |
| land | 128 | — | — | — | — |
| land | 512 | — | — | — | — |
| land | 2048 | — | — | — | — |
| land | 4608 | — | — | — | — |
| land | 8192 | — | — | — | — |
| land | 18432 | — | — | — | — |
| land | 41472 | — | — | — | — |
| land | 73728 | — | — | — | — |
| land | 165888 | — | — | — | — |
| land_no_vegetation | 128 | 318 | 1480 | — | — |
| land_no_vegetation | 512 | 83 | 1468 | — | — |
| land_no_vegetation | 2048 | 21 | 984 | — | — |
| land_no_vegetation | 4608 | 9.21 | 665 | — | — |
| land_no_vegetation | 8192 | 5.25 | 469 | — | — |
| land_no_vegetation | 18432 | 2.29 | 246 | — | — |
| land_no_vegetation | 41472 | 1.03 | 115 | — | — |
| land_no_vegetation | 73728 | 0.57 | 102 | — | — |
| land_no_vegetation | 165888 | 0.26 | 100 | — | — |
| soil_heat | 128 | 593 | 11470 | 85314 | 110298 |
| soil_heat | 512 | 152 | 9918 | 12579 | 111034 |
| soil_heat | 2048 | 38 | 6096 | 3151 | 94067 |
| soil_heat | 4608 | 17 | 16641 | 2346 | 62506 |
| soil_heat | 8192 | 9.51 | 3255 | 917 | 50495 |
| soil_heat | 18432 | 4.24 | 1637 | 415 | 21212 |
| soil_heat | 41472 | 1.87 | 760 | 150 | 11339 |
| soil_heat | 73728 | 1.06 | 320 | 93 | 7822 |
| soil_heat | 165888 | 0.47 | 595 | 42 | 4614 |

## Architecture: `cpu-x86`

Created for Terrarium.jl v0.1.4 on Sun, 09 Aug 2026 18:10:09 in `default` mode (1x time steps, 1 thread(s)).

### Machine details

```julia
julia> versioninfo()
Julia Version 1.12.2
Commit ca9b6662be4 (2025-11-20 16:25 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 128 × AMD EPYC 9554 64-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, znver4)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 128 virtual cores)
Environment:
  LD_LIBRARY_PATH = /usr/local/lib:/usr/local/lib:
```


### Model configurations, default resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land_no_vegetation | 3.75° | 4608 | 600 | 217 | 9.32 | 176 | 0.784 | 5.76 MiB | ok |
| soil_heat | 3.75° | 4608 | 600 | 217 | 17 | 96 | 1.44 | 3.85 MiB | ok |

### Land model, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 22.50° | 128 | — | — | — | — | — | — | failed: UndefVarError |
| land | 11.25° | 512 | — | — | — | — | — | — | failed: UndefVarError |
| land | 5.62° | 2048 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | 2.81° | 8192 | — | — | — | — | — | — | failed: UndefVarError |
| land | 1.88° | 18432 | — | — | — | — | — | — | failed: UndefVarError |
| land | 1.25° | 41472 | — | — | — | — | — | — | failed: UndefVarError |
| land | 0.94° | 73728 | — | — | — | — | — | — | failed: UndefVarError |
| land | 0.62° | 165888 | — | — | — | — | — | — | failed: UndefVarError |

### Land model without vegetation, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land_no_vegetation | 22.50° | 128 | 600 | 2000 | 318 | 5.16 | 0.744 | 171.16 KiB |
| land_no_vegetation | 11.25° | 512 | 600 | 1953 | 83 | 19.8 | 0.778 | 661.66 KiB |
| land_no_vegetation | 5.62° | 2048 | 600 | 488 | 21 | 79.3 | 0.774 | 2.56 MiB |
| land_no_vegetation | 3.75° | 4608 | 600 | 217 | 9.21 | 178 | 0.775 | 5.76 MiB |
| land_no_vegetation | 2.81° | 8192 | 600 | 122 | 5.25 | 313 | 0.785 | 10.23 MiB |
| land_no_vegetation | 1.88° | 18432 | 600 | 54 | 2.29 | 716 | 0.772 | 23.00 MiB |
| land_no_vegetation | 1.25° | 41472 | 600 | 24 | 1.03 | 1600 | 0.778 | 51.74 MiB |
| land_no_vegetation | 0.94° | 73728 | 600 | 20 | 0.57 | 2860 | 0.773 | 91.98 MiB |
| land_no_vegetation | 0.62° | 165888 | 600 | 20 | 0.26 | 6410 | 0.776 | 206.94 MiB |

### Soil heat conduction, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| soil_heat | 22.50° | 128 | 600 | 2000 | 593 | 2.77 | 1.39 | 114.63 KiB |
| soil_heat | 11.25° | 512 | 600 | 1953 | 152 | 10.8 | 1.43 | 443.13 KiB |
| soil_heat | 5.62° | 2048 | 600 | 488 | 38 | 43.3 | 1.42 | 1.72 MiB |
| soil_heat | 3.75° | 4608 | 600 | 217 | 17 | 97.3 | 1.42 | 3.85 MiB |
| soil_heat | 2.81° | 8192 | 600 | 122 | 9.51 | 173 | 1.42 | 6.85 MiB |
| soil_heat | 1.88° | 18432 | 600 | 54 | 4.24 | 388 | 1.43 | 15.40 MiB |
| soil_heat | 1.25° | 41472 | 600 | 24 | 1.87 | 877 | 1.42 | 34.65 MiB |
| soil_heat | 0.94° | 73728 | 600 | 20 | 1.06 | 1550 | 1.42 | 61.60 MiB |
| soil_heat | 0.62° | 165888 | 600 | 20 | 0.47 | 3490 | 1.43 | 138.59 MiB |

### Number of soil layers

| Configuration | Res | Columns | L | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | 10 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 20 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 30 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 60 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 100 | — | — | — | — | — | — | failed: UndefVarError |

### Number format, Float32 vs Float64

| Configuration | NF | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | Float32 | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | Float64 | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |

### Time stepper

| Configuration | Variant | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | ForwardEuler | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | Heun | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |

## Architecture: `gpu-nvidia`

Created for Terrarium.jl v0.1.4 on Sun, 09 Aug 2026 18:28:16 in `default` mode (1x time steps, 1 thread(s)).

### Machine details

```julia
julia> versioninfo()
Julia Version 1.12.2
Commit ca9b6662be4 (2025-11-20 16:25 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 128 × AMD EPYC 9554 64-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, znver4)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 128 virtual cores)
Environment:
  LD_LIBRARY_PATH = /usr/local/lib:/usr/local/lib:
```

```julia
julia> CUDA.versioninfo()
CUDA toolchain: 
- runtime 13.3.0, artifact installation
- driver 580.126.9 for 13.3
- compiler 13.3.33, artifact installation

CUDA libraries: 
- cuBLAS: 13.6.0
- cuSPARSE: 12.8.2
- cuSOLVER: 12.2.6
- cuFFT: 12.3.0
- cuRAND: 10.4.3
- CUPTI: 2026.2.1 (API 13.3.1)
- NVML: 13.0.0+580.126.9

Julia packages: 
- CUDACore: 6.2.1
- GPUArrays: 11.5.9
- GPUCompiler: 1.23.0
- KernelAbstractions: 0.9.42
- CUDA_Driver_jll: 13.3.0+1
- CUDA_Compiler_jll: 0.4.4+1
- CUDA_Runtime_jll: 0.23.0+1
- NVPTX_LLVM_Backend_jll: 22.1.7+1

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 78.157 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```


### Model configurations, default resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land_no_vegetation | 3.75° | 4608 | 600 | 217 | 678 | 2.42 | 57.1 | 5.76 MiB | ok |
| soil_heat | 3.75° | 4608 | 600 | 217 | 5387 | 0.305 | 453 | 3.85 MiB | ok |

### Land model, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 22.50° | 128 | — | — | — | — | — | — | failed: UndefVarError |
| land | 11.25° | 512 | — | — | — | — | — | — | failed: UndefVarError |
| land | 5.62° | 2048 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | 2.81° | 8192 | — | — | — | — | — | — | failed: UndefVarError |
| land | 1.88° | 18432 | — | — | — | — | — | — | failed: UndefVarError |
| land | 1.25° | 41472 | — | — | — | — | — | — | failed: UndefVarError |
| land | 0.94° | 73728 | — | — | — | — | — | — | failed: UndefVarError |
| land | 0.62° | 165888 | — | — | — | — | — | — | failed: UndefVarError |

### Land model without vegetation, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land_no_vegetation | 22.50° | 128 | 600 | 2000 | 1480 | 1.11 | 3.46 | 171.16 KiB |
| land_no_vegetation | 11.25° | 512 | 600 | 1953 | 1468 | 1.12 | 13.7 | 661.66 KiB |
| land_no_vegetation | 5.62° | 2048 | 600 | 488 | 984 | 1.67 | 36.8 | 2.56 MiB |
| land_no_vegetation | 3.75° | 4608 | 600 | 217 | 665 | 2.47 | 56 | 5.76 MiB |
| land_no_vegetation | 2.81° | 8192 | 600 | 122 | 469 | 3.5 | 70.1 | 10.23 MiB |
| land_no_vegetation | 1.88° | 18432 | 600 | 54 | 246 | 6.68 | 82.8 | 23.00 MiB |
| land_no_vegetation | 1.25° | 41472 | 600 | 24 | 115 | 14.3 | 86.8 | 51.74 MiB |
| land_no_vegetation | 0.94° | 73728 | 600 | 20 | 102 | 16.1 | 137 | 91.98 MiB |
| land_no_vegetation | 0.62° | 165888 | 600 | 20 | 100 | 16.5 | 303 | 206.94 MiB |

### Soil heat conduction, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| soil_heat | 22.50° | 128 | 600 | 2000 | 11470 | 0.143 | 26.8 | 114.63 KiB |
| soil_heat | 11.25° | 512 | 600 | 1953 | 9918 | 0.166 | 92.7 | 443.13 KiB |
| soil_heat | 5.62° | 2048 | 600 | 488 | 6096 | 0.269 | 228 | 1.72 MiB |
| soil_heat | 3.75° | 4608 | 600 | 217 | 16641 | 0.0987 | 1400 | 3.85 MiB |
| soil_heat | 2.81° | 8192 | 600 | 122 | 3255 | 0.505 | 487 | 6.85 MiB |
| soil_heat | 1.88° | 18432 | 600 | 54 | 1637 | 1 | 551 | 15.40 MiB |
| soil_heat | 1.25° | 41472 | 600 | 24 | 760 | 2.16 | 576 | 34.65 MiB |
| soil_heat | 0.94° | 73728 | 600 | 20 | 320 | 5.13 | 431 | 61.60 MiB |
| soil_heat | 0.62° | 165888 | 600 | 20 | 595 | 2.76 | 1800 | 138.59 MiB |

### Number of soil layers

| Configuration | Res | Columns | L | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | 10 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 20 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 30 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 60 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 100 | — | — | — | — | — | — | failed: UndefVarError |

### Number format, Float32 vs Float64

| Configuration | NF | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | Float32 | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | Float64 | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |

### Time stepper

| Configuration | Variant | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | ForwardEuler | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | Heun | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |

## Architecture: `reactant-cpu`

Created for Terrarium.jl v0.1.4 on Sun, 09 Aug 2026 18:53:09 in `default` mode (1.0x time steps, 1 thread(s)).

### Machine details

```julia
julia> versioninfo()
Julia Version 1.12.2
Commit ca9b6662be4 (2025-11-20 16:25 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 128 × AMD EPYC 9554 64-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, znver4)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 128 virtual cores)
Environment:
  LD_LIBRARY_PATH = /usr/local/lib:/usr/local/lib:
```

Reactant backend: `cpu` (selected with `Reactant.set_default_backend("cpu")`).


### Model configurations, default resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Compile | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | — | — | — | — | — | — | — | failed: UndefVarError |
| land_no_vegetation | 3.75° | 4608 | — | — | — | — | — | — | — | failed: ErrorException |
| soil_heat | 3.75° | 4608 | 600 | 217 | 1550 | 1.06 | 130 | 3.85 MiB | 31.5 s | ok |

### Land model, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 22.50° | 128 | — | — | — | — | — | — | failed: UndefVarError |
| land | 11.25° | 512 | — | — | — | — | — | — | failed: UndefVarError |
| land | 5.62° | 2048 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | 2.81° | 8192 | — | — | — | — | — | — | failed: UndefVarError |
| land | 1.88° | 18432 | — | — | — | — | — | — | failed: UndefVarError |
| land | 1.25° | 41472 | — | — | — | — | — | — | failed: UndefVarError |
| land | 0.94° | 73728 | — | — | — | — | — | — | failed: UndefVarError |
| land | 0.62° | 165888 | — | — | — | — | — | — | failed: UndefVarError |

### Land model without vegetation, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land_no_vegetation | 22.50° | 128 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 11.25° | 512 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 5.62° | 2048 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 3.75° | 4608 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 2.81° | 8192 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 1.88° | 18432 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 1.25° | 41472 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 0.94° | 73728 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 0.62° | 165888 | — | — | — | — | — | — | failed: ErrorException |

### Soil heat conduction, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| soil_heat | 22.50° | 128 | 600 | 2000 | 85314 | 0.0193 | 199 | 114.63 KiB | 14.6 s |
| soil_heat | 11.25° | 512 | 600 | 1953 | 12579 | 0.131 | 118 | 443.13 KiB | 14.6 s |
| soil_heat | 5.62° | 2048 | 600 | 488 | 3151 | 0.521 | 118 | 1.72 MiB | 14.6 s |
| soil_heat | 3.75° | 4608 | 600 | 217 | 2346 | 0.7 | 197 | 3.85 MiB | 1.1 s |
| soil_heat | 2.81° | 8192 | 600 | 122 | 917 | 1.79 | 137 | 6.85 MiB | 15.2 s |
| soil_heat | 1.88° | 18432 | 600 | 54 | 415 | 3.96 | 140 | 15.40 MiB | 15.7 s |
| soil_heat | 1.25° | 41472 | 600 | 24 | 150 | 10.9 | 114 | 34.65 MiB | 16.8 s |
| soil_heat | 0.94° | 73728 | 600 | 20 | 93 | 17.7 | 125 | 61.60 MiB | 18.5 s |
| soil_heat | 0.62° | 165888 | 600 | 20 | 42 | 39.1 | 127 | 138.59 MiB | 23.0 s |

### Number of soil layers

| Configuration | Res | Columns | L | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | 10 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 20 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 30 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 60 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 100 | — | — | — | — | — | — | failed: UndefVarError |

### Number format, Float32 vs Float64

| Configuration | NF | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | Float32 | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | Float64 | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |

### Time stepper

| Configuration | Variant | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | ForwardEuler | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | Heun | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |

## Architecture: `reactant-gpu`

Created for Terrarium.jl v0.1.4 on Sun, 09 Aug 2026 18:38:21 in `default` mode (1x time steps, 1 thread(s)).

### Machine details

```julia
julia> versioninfo()
Julia Version 1.12.2
Commit ca9b6662be4 (2025-11-20 16:25 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 128 × AMD EPYC 9554 64-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, znver4)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 128 virtual cores)
Environment:
  LD_LIBRARY_PATH = /usr/local/lib:/usr/local/lib:
```

```julia
julia> CUDA.versioninfo()
CUDA toolchain: 
- runtime 13.3.0, artifact installation
- driver 580.126.9 for 13.3
- compiler 13.3.33, artifact installation

CUDA libraries: 
- cuBLAS: 13.6.0
- cuSPARSE: 12.8.2
- cuSOLVER: 12.2.6
- cuFFT: 12.3.0
- cuRAND: 10.4.3
- CUPTI: 2026.2.1 (API 13.3.1)
- NVML: 13.0.0+580.126.9

Julia packages: 
- CUDACore: 6.2.1
- GPUArrays: 11.5.9
- GPUCompiler: 1.23.0
- KernelAbstractions: 0.9.42
- CUDA_Driver_jll: 13.3.0+1
- CUDA_Compiler_jll: 0.4.4+1
- CUDA_Runtime_jll: 0.23.0+1
- NVPTX_LLVM_Backend_jll: 22.1.7+1

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 19.239 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```

Reactant backend: `gpu` (selected with `Reactant.set_default_backend("gpu")`).


### Model configurations, default resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Compile | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | — | — | — | — | — | — | — | failed: UndefVarError |
| land_no_vegetation | 3.75° | 4608 | — | — | — | — | — | — | — | failed: ErrorException |
| soil_heat | 3.75° | 4608 | 600 | 217 | 65324 | 0.0251 | 5500 | 3.85 MiB | 31.7 s | ok |

### Land model, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 22.50° | 128 | — | — | — | — | — | — | failed: UndefVarError |
| land | 11.25° | 512 | — | — | — | — | — | — | failed: UndefVarError |
| land | 5.62° | 2048 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | 2.81° | 8192 | — | — | — | — | — | — | failed: UndefVarError |
| land | 1.88° | 18432 | — | — | — | — | — | — | failed: UndefVarError |
| land | 1.25° | 41472 | — | — | — | — | — | — | failed: UndefVarError |
| land | 0.94° | 73728 | — | — | — | — | — | — | failed: UndefVarError |
| land | 0.62° | 165888 | — | — | — | — | — | — | failed: UndefVarError |

### Land model without vegetation, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land_no_vegetation | 22.50° | 128 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 11.25° | 512 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 5.62° | 2048 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 3.75° | 4608 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 2.81° | 8192 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 1.88° | 18432 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 1.25° | 41472 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 0.94° | 73728 | — | — | — | — | — | — | failed: ErrorException |
| land_no_vegetation | 0.62° | 165888 | — | — | — | — | — | — | failed: ErrorException |

### Soil heat conduction, horizontal resolution

| Configuration | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| soil_heat | 22.50° | 128 | 600 | 2000 | 110298 | 0.0149 | 258 | 114.63 KiB | 14.7 s |
| soil_heat | 11.25° | 512 | 600 | 1953 | 111034 | 0.0148 | 1040 | 443.13 KiB | 14.6 s |
| soil_heat | 5.62° | 2048 | 600 | 488 | 94067 | 0.0175 | 3520 | 1.72 MiB | 14.8 s |
| soil_heat | 3.75° | 4608 | 600 | 217 | 62506 | 0.0263 | 5260 | 3.85 MiB | 1.1 s |
| soil_heat | 2.81° | 8192 | 600 | 122 | 50495 | 0.0325 | 7550 | 6.85 MiB | 15.8 s |
| soil_heat | 1.88° | 18432 | 600 | 54 | 21212 | 0.0774 | 7140 | 15.40 MiB | 15.7 s |
| soil_heat | 1.25° | 41472 | 600 | 24 | 11339 | 0.145 | 8590 | 34.65 MiB | 16.7 s |
| soil_heat | 0.94° | 73728 | 600 | 20 | 7822 | 0.21 | 10500 | 61.60 MiB | 18.3 s |
| soil_heat | 0.62° | 165888 | 600 | 20 | 4614 | 0.356 | 14000 | 138.59 MiB | 23.0 s |

### Number of soil layers

| Configuration | Res | Columns | L | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | 3.75° | 4608 | 10 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 20 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 30 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 60 | — | — | — | — | — | — | failed: UndefVarError |
| land | 3.75° | 4608 | 100 | — | — | — | — | — | — | failed: UndefVarError |

### Number format, Float32 vs Float64

| Configuration | NF | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | Float32 | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | Float64 | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |

### Time stepper

| Configuration | Variant | Res | Columns | Δt | Steps | SYPD | ms/step | Mcell-steps/s | Memory | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| land | ForwardEuler | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |
| land | Heun | 3.75° | 4608 | — | — | — | — | — | — | failed: UndefVarError |

