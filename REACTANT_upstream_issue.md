# Upstream issue draft: Oceananigans — `RectilinearGrid` with explicit array vertical coordinates fails Reactant kernel tracing

**Where to file:** Oceananigans.jl. Only needed for **non-uniform** vertical grids; regular-range
`z` works today, so this is not urgent for uniform-spacing use.

## Summary

A `HydrostaticFreeSurfaceModel` (or any halo-filling op) on a `RectilinearGrid` compiles and runs
under Reactant when the vertical coordinate is a regular **range** (`z=(zmin,zmax)` →
`StepRangeLen`), but **fails during compilation** when the vertical coordinate is an explicit
**array/`Vector`** (e.g. stretched/exponential grids). The grid's coordinate arrays are not
converted to `CuTracedArray` during kernel adaptation, and Reactant's
`_check_no_traced_in_kernel_arg` rejects the grid kernel argument.

The Reactant extension already handles array coordinates for `LatitudeLongitudeGrid` (via a
dedicated `on_architecture`/`_to_reactant` path) but not for `RectilinearGrid`; a source comment
notes RG "retains StepRangeLen coordinates which sidestep the CuTracedArray VX type constraint",
i.e. array-coordinate RG is effectively unsupported.

## Environment

- Julia 1.12.6 (aarch64, macOS), CPU backend (`Reactant.set_default_backend("cpu")`, `CUDA` loaded)
- Oceananigans **0.110.7** (latest), Reactant **0.2.270**, KernelAbstractions 0.9.42

## Reproduction (A/B — same model, only `z` differs)

```julia
using Oceananigans, Reactant, CUDA
using Oceananigans.Architectures: ReactantState
using Reactant: @trace
arch = ReactantState()

run4!(m) = begin
    f(m,Δt,N) = (@trace track_numbers=false for _ in 1:N; time_step!(m,Δt); end; nothing)
    g = @compile raise=true raise_first=true sync=true f(m, 60.0, 4); g(m, 60.0, 4)
end

# A) regular range z  -> COMPILES AND RUNS
gA = RectilinearGrid(arch; size=16, z=(-200, 0), topology=(Flat, Flat, Bounded))
run4!(HydrostaticFreeSurfaceModel(gA; buoyancy=nothing, tracers=:T))   # OK

# B) explicit z vector -> FAILS during compilation
zv = collect(range(-200f0, 0f0, length=17))
gB = RectilinearGrid(arch, Float32; size=16, z=zv, topology=(Flat, Flat, Bounded))
run4!(HydrostaticFreeSurfaceModel(gB; buoyancy=nothing, tracers=:T))   # error (see below)
```

## Observed (case B)

```
GPU kernel argument of type RectilinearGrid{TracedRNumber{Float32}, …,
  StaticVerticalDiscretization{OffsetVector{TracedRNumber{Float32}, TracedRArray{Float32,1}}, …}, …}
All TracedRNumber/TracedRArray must have been replaced by their CuTracedRNumber/CuTracedArray …
```

## Expected

`RectilinearGrid` with array vertical coordinates should adapt its coordinate arrays to
`CuTracedArray` for kernel launch, as `LatitudeLongitudeGrid` already does — so stretched vertical
grids trace under Reactant.

## Notes

- Downstream (Terrarium.jl) works around this for **uniform** spacing by passing a range `z`; only
  stretched/exponential/prescribed vertical grids need this upstream fix.
- Likely fix location: the `on_architecture(::ReactantState, ::RectilinearGrid)` /
  `StaticVerticalDiscretization` adapt path in `OceananigansReactantExt`, mirroring the
  `LatitudeLongitudeGrid` treatment.
