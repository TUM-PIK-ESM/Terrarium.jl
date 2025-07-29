"""
Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Fixed timestep size in seconds"
    dt::NF = 300.0
end

get_dt(euler::ForwardEuler) = euler.dt

is_adaptive(euler::ForwardEuler) = false
