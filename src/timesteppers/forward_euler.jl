"""
Simple forward Euler time stepping scheme.
"""
@kwdef struct ForwardEuler{NF} <: AbstractTimeStepper{NF}
    "Fixed timestep size in seconds"
    dt::NF = 300.0
end

struct ForwardEulerCache{NF} <: AbstractTimeStepperCache{NF} end

get_dt(euler::ForwardEuler) = euler.dt

is_adaptive(euler::ForwardEuler) = false

@inline step(::ForwardEuler, progvar, tendency, dt) = progvar + dt*tendency

initialize(::AbstractModel{NF}, ::ForwardEuler) where {NF} = ForwardEulerCache{NF}()

# function timestep!(state, model::AbstractModel, euler::ForwardEuler, dt=get_dt(euler))
#     # TODO: implement timestep! generically
#     error("not implemented")
# end
