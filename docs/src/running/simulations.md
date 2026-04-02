# Running simulations

```@meta
CurrentModule = Terrarium
```

AS we have already seen, running a simulation in Terrarium involves a minimum of three steps:

1. **Build a model** — construct an `AbstractModel` subtype with a `grid`, processes, and constants.
2. **Initialize an integrator** — call `initialize` to allocate state variables and prepare the time stepper.
3. **Run** — advance the simulation via `run!`, or step it manually with `timestep!`.

So far, however, these steps only allow us to run transient, open-ended simulations. For our model to be useful,
we need to be able to define stopping criteria, diagnostics, and outputs. Fortunately, since Terrarium is built
on top of `Oceananigans`, we can make use of one of its core utilities that can help us define simulations:

```@docs
Simulation
```

## Setting up a `Simulation`

### Building a model

All Terrarium simulations begin with a model. We will use again here [`SoilModel`](@ref) as an example:
```@example simulation
using Terrariumstepper

grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 10))
model = SoilModel(grid)
```

As discussed in the section on [initialization](@ref "Initialization"), [`initialize`](@ref) allocates
all state variables, runs the model [`initialize!`](@ref) routine, and wraps everything in a [`ModelIntegrator`](@ref):
```@example
integrator = initialize(model, Heun(eltype(grid));
    initializers = (temperature = 5.0,),     # set individual field initial values
)
```

### Choosing a time stepper

Terrarium currently provides two choices of explicit time steppers: [`ForwardEuler`](@ref) and [`Heun`](@ref).

```@docs; canonical
ForwardEuler
Heun
```

Both are constructed with a default timestep `Δt` in seconds. `Δt` can be overridden at run
time by specifying it when calling `run!` or `timestep!`.

### Transient simulations

The [`ModelIntegrator`](@ref) is already runnable without an `Oceananigans.Simulation` wrapper:

```@example simulation
# Advance by a fixed number of steps
run!(integrator; steps = 100)

# Advance for a given wall-clock period
using Dates: Day

run!(integrator; period = Day(30), Δt = 3600.0)

# Advance a single step
timestep!(integrator)
```

Recall that state variable `Field`s can be accessed via `integrator.state`:

```@example simulation
integrator.state.temperature   # current temperature Field
```

### Setting up a `Simulation`

Since Terrarium's [`ModelIntegrator`](@ref) implements the Oceananigans model interface, we can
wrap the integrator in a [`Simulation`]:

```@example simulation
using Oceananigans.Units: days

sim = Simulation(integrator; stop_time = 30days, Δt = 300.0)
run!(sim)
```

Note that `stop_time` is usually expressed in the same units as the `clock` — seconds by default.
However, `Oceananigans.Units` provides convenient conversion factors (e.g. `days`, `hours`).

!!! note
    Always call `Terrarium.initialize!(integrator)` **before** attaching output writers and
    calling `run!(sim)` when re-running a simulation; otherwise, the clock and state may already
    be at the end of the simulation and output writers will not trigger.

## Output writers

We can make use of `Oceananigans.OutputWriters` to flexibly save output to disk in a variety of formats.

The simplest choice output writer is `JLD2Writer` which saves selected `Field`s to a `.jld2` file at a specified schedule:

```@example simulation
using Oceananigans: JLD2Writer, TimeInterval
using Oceananigans.Units: hours

sim.output_writers[:soil] = JLD2Writer(
    integrator,
    (temperature = integrator.state.temperature,
     saturation  = integrator.state.saturation_water_ice);
    filename = "output/soil_run.jld2",
    overwrite_existing = true,
    schedule = TimeInterval(6hours),
)
```

The saved file can be read back as a `FieldTimeSeries` for post-processing:

```@example simulation
using Oceananigans: FieldTimeSeries

temperature_ts = FieldTimeSeries("output/soil_run.jld2", "temperature")
temperature_ts[end]   # Field at the last saved time
```

Output writers accept any Oceananigans schedule as the `schedule` keyword:

| Schedule type | Description |
|---------------|-------------|
| `TimeInterval(Δt)` | Write every `Δt` seconds of simulation time |
| `IterationInterval(n)` | Write every `n` timesteps |
| `AveragedTimeInterval(Δt)` | Write time-averaged outputs over windows of `Δt` seconds |

These are all re-exported from Oceananigans and require no additional imports beyond
`using Oceananigans`.

Multiple output writers can be added to the same simulation, e.g. to save different
variables at different frequencies:

```julia
sim.output_writers[:fast]   = JLD2Writer(integrator, (temperature = ...,); schedule = TimeInterval(1hours),  ...)
sim.output_writers[:daily]  = JLD2Writer(integrator, (pressure_head = ...,); schedule = AveragedTimeInterval(1days), ...)
```

## Callbacks

Terrarium simulations inherit the full Oceananigans `Callback` machinery.
A callback is a function that is called by `Simulation` at a given schedule during `run!`:

```julia
using Oceananigans: Callback, IterationInterval

function print_progress(sim)
    t = sim.model.clock.time
    iter = sim.model.clock.iteration
    @info "Iteration $iter, time $t s"
end

sim.callbacks[:progress] = Callback(print_progress, IterationInterval(100))
```

The callback function receives the `Simulation` object, giving it access to the full
`integrator` via `sim.model`, and to the current `state` via `sim.model.state`.
