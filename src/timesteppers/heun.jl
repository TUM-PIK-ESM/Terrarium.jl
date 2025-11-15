"""
    $TYPEDEF

Simple forward 2nd order Heun / improved Euler time stepping scheme.
"""
@kwdef struct Heun{NF,S} <: AbstractTimeStepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0

    "Stage (Cache) for intermediate results"
    stage::S
end

"""
    Heun(state::StateVariables, Δt=300)

Create a `Heun` timestepper for the given state variables.
"""
Heun(state::StateVariables; kwargs...) = Heun{eltype(state), typeof(state)}(; stage=deepcopy(state), kwargs...)

default_dt(heun::Heun) = heun.Δt

is_adaptive(heun::Heun) = false

function average_tendencies!(
    state::StateVariables{NF, prognames, tendnames},
    stage::StateVariables{NF, prognames, tendnames}
) where {NF, prognames, tendnames}
    for tendname in tendnames
        state.tendencies[tendname] .= (state.tendencies[tendname] + stage.tendencies[tendname]) / 2
    end
end

function timestep!(state, model::AbstractModel, timestepper::Heun, Δt = default_dt(timestepper))
    # Copy current state to stage
    stage = timestepper.stage
    copyto!(stage, state)

    # Compute stage
    # TODO: this is currently incorrect because inputs are not updated
    explicit_step!(stage, get_grid(model), timestepper, Δt)
    reset_tendencies!(stage)
    compute_auxiliary!(stage, model) 
    compute_tendencies!(stage, model) 

    # Final improved Euler step call that steps `state` forward but averages `state.tendencies`
    average_tendencies!(state, stage)
    explicit_step!(state, get_grid(model), timestepper, Δt) 
    
    # Apply inverse closure relations
    for closure in state.closures
        invclosure!(state, model, closure)
    end

    # Update clock and return
    tick!(state.clock, Δt)
    return nothing
end
