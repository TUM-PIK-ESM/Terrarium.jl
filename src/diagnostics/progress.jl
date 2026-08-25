"""
    $TYPEDEF

Simple callable `struct` wrapper around `ProgressMeter.Progress` intended to be
used in an Oceananigans `Callback`.

```julia
simulation.callbacks[:progress] = Callback(ProgressReporter(nsteps), schedule)
```

The callback can be created directly via [`ProgressCallback`](@ref),

```julia
simulation.callbacks[:progress] = ProgressCallback(nsteps)
```
"""
struct ProgressReporter{Progress <: ProgressMeter.Progress}
    progress::Progress
end

ProgressReporter(n::Int = 100; kwargs...) = ProgressReporter(ProgressMeter.Progress(n; kwargs...))
function ProgressReporter(simulation::Simulation; kwargs...)
    desc = "Running simulation (Δt = $(simulation.Δt))"
    if isfinite(simulation.stop_iteration)
        return ProgressReporter(Int(simulation.stop_iteration); desc, kwargs...)
    else
        return ProgressReporter(; kwargs...)
    end
end

function (callback::ProgressReporter)(simulation)
    desc = "Running simulation (Δt = $(simulation.Δt))"
    if isfinite(simulation.stop_iteration)
        return ProgressMeter.update!(callback.progress, simulation.model.clock.iteration; desc)
    else
        val = Int(ceil(simulation.model.clock.time / simulation.stop_time * 100))
        return ProgressMeter.update!(callback.progress, val; desc)
    end
end

"""
Creates a `Callback` for `ProgressReporter` with schedule `IterationInterval(1)`.
"""
ProgressCallback(args...; kwargs...) = Callback(ProgressReporter(args...; kwargs...), IterationInterval(1))
