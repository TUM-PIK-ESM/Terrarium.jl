"""
Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Fixed timestep size in seconds"
    Δt::NF = 300.0
end

ForwardEuler(::Type{NF}; kwargs...) where {NF} = ForwardEuler{NF}(; kwargs...)

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
