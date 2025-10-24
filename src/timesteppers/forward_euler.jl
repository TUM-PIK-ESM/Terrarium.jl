"""
    $TYPEDEF

Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0
end

"""
    ForwardEuler(state::StateVariables; kwargs...)

Create a `ForwardEuler` timestepper with the given numeric format `NF`.
"""
ForwardEuler(state::StateVariables; kwargs...) = ForwardEuler{eltype(state)}(; kwargs...)

default_dt(euler::ForwardEuler) = euler.Δt

is_adaptive(euler::ForwardEuler) = false

function timestep!(state, model::AbstractModel, timestepper::ForwardEuler, Δt = default_dt(timestepper))
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    # Euler step
    explicit_step!(state, get_grid(model), timestepper, Δt)
    # Apply inverse closure relations
    for closure in state.closures
        invclosure!(state, model, closure)
    end
    # Update clock
    tick!(state.clock, Δt)
    return nothing
end
