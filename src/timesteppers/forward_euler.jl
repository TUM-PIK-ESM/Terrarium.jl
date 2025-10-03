"""
Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Fixed timestep size in seconds"
    dt::NF = 300.0
end

get_dt(euler::ForwardEuler) = euler.dt

is_adaptive(euler::ForwardEuler) = false

initialize(::AbstractModel{NF}, ::ForwardEuler) where {NF} = nothing

function timestep!(state, model::AbstractModel, timestepper::AbstractTimeStepper, dt = get_dt(timestepper))
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    # perform Euler step
    explicit_step!(state, model, timestepper, dt)
end
