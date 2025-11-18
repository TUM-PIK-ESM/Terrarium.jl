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
    state::StateVariables{NF, prognames, tendnames},
    stage::StateVariables{NF, prognames, tendnames}
) where {NF, prognames, tendnames}
    for tendname in tendnames
        state.tendencies[tendname] .= (state.tendencies[tendname] + stage.tendencies[tendname]) / 2
    end
end

function timestep!(state, timestepper::Heun, model::AbstractModel, inputs::InputSources, Δt = default_dt(timestepper))
    @assert is_initialized(timestepper)

    grid = get_grid(model)

    # Copy current state to stage
    stage = timestepper.stage
    copyto!(stage, state)

    # Compute stage
    explicit_step!(stage, grid, timestepper, Δt)
    # Apply closure relations
    closure!(stage, model)
    # Update clock
    tick!(stage.clock, Δt)
    update_state!(stage, model, inputs, compute_tendencies = true)

    # Final improved Euler step call that steps `state` forward but averages `state.tendencies`
    average_tendencies!(state, stage)
    explicit_step!(state, grid, timestepper, Δt) 
    # Apply closure relations
    closure!(state, model)
    # Update clock
    tick!(state.clock, Δt)
    return nothing
end
