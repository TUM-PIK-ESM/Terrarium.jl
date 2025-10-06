"""
Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Fixed timestep size in seconds"
    dt::NF = 300.0
end

default_dt(euler::ForwardEuler) = euler.dt

is_adaptive(euler::ForwardEuler) = false

function timestep!(state, model::AbstractModel, timestepper::ForwardEuler, dt = default_dt(timestepper))
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    # Euler step
    explicit_step!(state, get_grid(model), timestepper, dt)
    # Apply inverse closure relations
    for closure in state.closures
        invclosure!(state, model, closure)
    end
    # Update clock
    tick!(state.clock, dt)
    return nothing
end
