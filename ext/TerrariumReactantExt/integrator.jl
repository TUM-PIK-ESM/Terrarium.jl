# Compiled time stepping.
#
# `run!` on a device integrator compiles the whole stepping loop into a single StableHLO
# program (via `@trace for`) and runs it. The compiled function can be reused across runs by
# passing it back as the `compiled_run!` keyword. The compiled core calls the timestepper-level
# `timestep!` (`timestep!(integrator, timestepper, Δt)`), which is NOT overridden here — this is
# what avoids recursively re-entering the compiling wrapper.

const ReactantIntegrator = ModelIntegrator{<:Any, <:RARCH}

# One full step + auxiliary finalization; the unit stepped by the traced loop.
function step_core!(integrator, Δt)
    timestepper = get_timestepper(integrator.model)
    Terrarium.timestep!(integrator, timestepper, Terrarium.convert_dt(Δt))
    Terrarium.compute_auxiliary!(integrator.state, integrator.model)
    return nothing
end

# Advance the integrator by `Nt` steps of size `Δt`. `Reactant.@trace` compiles the loop to a
# single `stablehlo.while` (rather than unrolling `Nt` copies of the step into the trace). This
# is the function compiled by `run!` and mirrors Oceananigans' `run_timesteps!` test helper.
function run_timesteps!(integrator, Δt, Nt)
    Reactant.@trace track_numbers = false for _ in 1:Nt
        step_core!(integrator, Δt)
    end
    return nothing
end

"""
    run!(integrator; steps, period, Δt, compiled_run! = nothing)

Run a `ReactantState` integrator. If `compiled_run!` is `nothing`, the stepping loop is compiled
here with Reactant (`raise=true` lowers the KernelAbstractions kernels — halo fills,
tendency/closure kernels — so their traced grid/array arguments are adapted correctly).
The compiled function is returned-through by being reusable: pass it back as `compiled_run!` on
subsequent runs with the same `Δt` and step count to skip recompilation.
"""
function Oceananigans.Simulations.run!(
        integrator::ReactantIntegrator;
        steps::Union{Int, Nothing} = nothing,
        period::Union{Terrarium.Period, Nothing} = nothing,
        Δt = Terrarium.default_dt(get_timestepper(integrator.model)),
        compiled_run! = nothing,
    )
    Δt = Terrarium.convert_dt(Δt)
    nsteps = Terrarium.get_steps(steps, period, Δt)
    if isnothing(compiled_run!)
        @info "Reactant: compiling run_timesteps! (Δt=$Δt, steps=$nsteps)"
        compiled_run! = @compile raise = true raise_first = true sync = true run_timesteps!(integrator, Δt, nsteps)
    end
    compiled_run!(integrator, Δt, nsteps)
    return integrator
end

# Single step: the same compiled path with one iteration (no persistent compilation cache — a
# repeated single-step loop should use `run!` instead, which compiles the loop once).
function Terrarium.timestep!(integrator::ReactantIntegrator, Δt; finalize = true)
    run!(integrator; steps = 1, Δt)
    return nothing
end
