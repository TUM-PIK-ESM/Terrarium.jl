# Compiled time stepping.
#
# `timestep!`/`run!` on a device integrator compile a single step with Reactant (once, cached) and
# drive it from a host loop. The compiled core calls the timestepper-level `timestep!`
# (`timestep!(integrator, timestepper, Δt)`), which is NOT overridden here — this is what avoids
# recursively re-entering the compiling wrapper.

const ReactantIntegrator = ModelIntegrator{<:Any, <:RARCH}

# Compiled steps keyed by (integrator identity, Δt). Δt is a trace constant, so a new Δt recompiles.
const COMPILED_STEPS = Dict{Any, Any}()

# One full step + auxiliary finalization; this is what gets traced/compiled.
function step_core!(integrator, Δt)
    timestepper = get_timestepper(integrator.model)
    Terrarium.timestep!(integrator, timestepper, Terrarium.convert_dt(Δt))
    Terrarium.compute_auxiliary!(integrator.state, integrator.model)
    return nothing
end

function compiled_step(integrator::ReactantIntegrator, Δt)
    return get!(COMPILED_STEPS, (objectid(integrator), Δt)) do
        @info "Reactant: compiling timestep! (Δt=$Δt)"
        # `raise=true` lowers the KernelAbstractions kernels (halo fills, tendency/closure kernels)
        # so their traced grid/array arguments are adapted correctly; without it the launch trips
        # Reactant's `_check_no_traced_in_kernel_arg`. Mirrors Oceananigans' own Reactant tests.
        @compile raise = true raise_first = true sync = true step_core!(integrator, Δt)
    end
end

function Terrarium.timestep!(integrator::ReactantIntegrator, Δt; finalize = true)
    r_step! = compiled_step(integrator, Δt)
    r_step!(integrator, Δt)
    return nothing
end

function Oceananigans.Simulations.run!(
        integrator::ReactantIntegrator;
        steps::Union{Int, Nothing} = nothing,
        period::Union{Terrarium.Period, Nothing} = nothing,
        Δt = Terrarium.default_dt(get_timestepper(integrator.model)),
    )
    Δt = Terrarium.convert_dt(Δt)
    nsteps = Terrarium.get_steps(steps, period, Δt)
    r_step! = compiled_step(integrator, Δt)
    for _ in 1:nsteps
        r_step!(integrator, Δt)
    end
    return integrator
end
