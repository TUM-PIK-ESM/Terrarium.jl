"""
Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Fixed timestep size in seconds"
    dt::NF = 300.0
end

get_dt(euler::ForwardEuler) = euler.dt

is_adaptive(euler::ForwardEuler) = false

function timestep!(state, model::AbstractModel, euler::ForwardEuler, dt=get_dt(euler))
    # TODO: implement timestep! generically
    error("not implemented")
end
