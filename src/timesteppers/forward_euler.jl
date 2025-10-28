"""
    $TYPEDEF

Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0
end

"""
    ForwardEuler(::Type{NF}; kwargs...)

Create a `ForwardEuler` timestepper with the given numeric format `NF`.
"""
ForwardEuler(::Type{NF}; kwargs...) where {NF} = ForwardEuler{NF}(; kwargs...)

default_dt(euler::ForwardEuler) = euler.Δt

is_adaptive(euler::ForwardEuler) = false

function timestep!(state, model::AbstractModel, timestepper::ForwardEuler, Δt = default_dt(timestepper))
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    do_timestep!(state, model, timestepper, Δt)
    return nothing
end

function do_timestep!(state, model::AbstractModel, timestepper::ForwardEuler, Δt)
    # Euler step
    explicit_step!(state, get_grid(model), timestepper, Δt)
    # Apply inverse closure relations
    invclosure!(state, model)
    # Update clock
    tick!(state.clock, Δt)
end
