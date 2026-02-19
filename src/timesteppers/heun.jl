"""
    $TYPEDEF

Simple forward 2nd order Heun / improved Euler time stepping scheme.
"""
@kwdef struct Heun{NF, Stage} <: AbstractTimeStepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0

    "Stage (Cache) for intermediate results"
    stage::Stage = nothing
end

default_dt(heun::Heun) = heun.Δt

is_adaptive(heun::Heun) = false

is_initialized(heun::Heun) = !isnothing(heun.stage)

function initialize(timestepper::Heun, ::AbstractModel, state)
    return setproperties(timestepper, stage = deepcopy(state))
end

function average_tendencies!(
        state::StateVariables{NF, prognames},
        stage::StateVariables{NF, prognames}
    ) where {NF, prognames}
    for name in prognames
        state.tendencies[name] .= (state.tendencies[name] + stage.tendencies[name]) / 2
    end
    return
end

function timestep!(integrator::ModelIntegrator, timestepper::Heun, Δt = default_dt(timestepper))
    @assert is_initialized(timestepper)

    (; model, state, inputs) = integrator
    grid = get_grid(model)

    # Copy current state to stage
    stage = timestepper.stage
    copyto!(stage, state)

    # Compute stage
    update_state!(stage, model, inputs, compute_tendencies = true)
    explicit_step!(stage, grid, timestepper, Δt)
    # Call timestep! on model
    timestep!(stage, model, timestepper, Δt)
    # Apply closure relations
    closure!(stage, model)
    # Update clock
    tick!(stage.clock, Δt)
    # Recompute tendencies after timestep
    update_state!(stage, model, inputs, compute_tendencies = true)

    # Final improved Euler step call that steps `state` forward but averages `state.tendencies`
    average_tendencies!(state, stage)
    explicit_step!(state, grid, timestepper, Δt)
    # Call timestep! on model
    timestep!(state, model, timestepper, Δt)
    # Apply closure relations
    closure!(state, model)
    # Update clock
    tick!(state.clock, Δt)
    return nothing
end
