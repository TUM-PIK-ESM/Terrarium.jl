"""
    $TYPEDEF

Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0
end

ForwardEuler(::Type{NF}; kwargs...) where {NF} = ForwardEuler{NF}(; kwargs...)

timestepping(::ForwardEuler) = Explicit()

default_dt(euler::ForwardEuler) = euler.Δt

is_adaptive(euler::ForwardEuler) = false

# Step only the prognostic variables named in `names` with the forward Euler scheme.
function timestep!(integrator::ModelIntegrator, timestepper::ForwardEuler, Δt, names::Tuple)
    # Compute auxiliaries and tendencies
    update_state!(integrator, compute_tendencies = true)
    # Euler step for the assigned prognostic variables
    explicit_step!(integrator.state, get_grid(integrator.model), timestepper, Δt, names)
    # Call timestep! on model
    timestep!(integrator.state, integrator.model, timestepper, Δt)
    # Apply closure relations
    closure!(integrator.state, integrator.model)
    return nothing
end
