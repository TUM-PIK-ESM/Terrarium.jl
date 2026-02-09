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

function timestep!(state, timestepper::ForwardEuler, model::AbstractModel, inputs::InputSources, Δt)
    # Euler step
    explicit_step!(state, get_grid(model), timestepper, Δt)
    # Apply closure relations
    closure!(state, model)
    # Update clock
    return tick!(state.clock, Δt)
end
