# Accelerating simulations with Reactant

```@meta
CurrentModule = Terrarium
```

Terrarium models can be compiled and executed with
[Reactant.jl](https://github.com/EnzymeAD/Reactant.jl), which traces the model time step into
MLIR/StableHLO and compiles it with XLA for high-performance execution on CPUs, GPUs, and TPUs.
The only user-facing change is the architecture passed to the grid: with `ReactantState()`,
initialization runs on the CPU and is transferred to the device, and the first call to
`timestep!` or `run!` compiles the time step transparently (subsequent steps reuse the compiled
program).

!!! warning "Experimental"
    Reactant support is experimental and currently validated for `SoilModel` heat conduction
    on `ColumnGrid` and `ColumnRingGrid` with **uniform vertical spacing**. Correctness against
    the CPU implementation is tested continuously in `test/reactant/`.

## Example

The following mirrors the [soil heat conduction example](@ref soil_heat_column), changing only
the architecture. Note that this code is not executed during the documentation build since
Reactant compilation takes several minutes.

```julia
using Terrarium
using Reactant, CUDA  # CUDA is required by Reactant's kernel integration, even on CPU

# ReactantState() instead of CPU() or GPU() — the only change!
grid = ColumnGrid(ReactantState(), Float32, UniformSpacing(Δz = 0.2f0, N = 10))
model = SoilModel(grid)

boundary_conditions = PrescribedSurfaceTemperature(:T_ub, 1.0f0)
initializers = (temperature = (x, z) -> -1.0f0 - 0.05f0 * z,)
integrator = initialize(model; boundary_conditions, initializers)

# The first step compiles the model with XLA (takes a moment); further steps are fast.
run!(integrator, steps = 144, Δt = 600.0f0)

# Materialize device results on the CPU for analysis/plotting
T = Array(interior(integrator.state.temperature))
```

The target device is selected by Reactant, e.g. `Reactant.set_default_backend("gpu")` before
constructing the model.

In the example folder you can also find examples that demonstrate how to use the Reactant model to take derivatives of the model and integratate and train neural networks. 