"""
    $TYPEDEF

Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0
end

ForwardEuler(::Type{NF}; kwargs...) where {NF} = ForwardEuler{NF}(; kwargs...)

default_dt(euler::ForwardEuler) = euler.Δt

is_adaptive(euler::ForwardEuler) = false

is_initialized(euler::ForwardEuler) = true

function timestep!(integrator::ModelIntegrator, timestepper::ForwardEuler, Δt)
    # Compute auxiliaries and tendencies
    update_state!(integrator, compute_tendencies = true)
    # Euler step
    explicit_step!(integrator.state, get_grid(integrator.model), timestepper, Δt)
    # Call timestep! on model
    timestep!(integrator.state, integrator.model, timestepper, Δt)
    # Apply closure relations
    closure!(integrator.state, integrator.model)
    # Update clock
    tick!(integrator.state.clock, Δt)
    return nothing
end
