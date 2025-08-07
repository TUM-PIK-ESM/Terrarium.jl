abstract type AbstractSimulation end

struct Simulation{
    Model<:AbstractModel,
    TimeStepper<:AbstractTimeStepper,
    StateVars<:AbstractStateVariables,
} <: AbstractSimulation
    "The type of model used for the simulation."
    model::Model

    "The time stepping scheme used by the simulation."
    time_stepping::TimeStepper

    "Collection of all state variables defined on the simulation `model`."
    state::StateVars
end

function initialize(model::AbstractModel, time_stepping=get_time_stepping(model); clock::Clock=Clock(time=0.0))
    state = StateVariables(model, clock)
    sim = Simulation(model, time_stepping, state)
    initialize!(sim)
    return sim
end

function initialize!(sim::Simulation)
    # TODO: reset other variables too?
    reset_tendencies!(sim.state)
    initialize!(sim.state, sim.model)
    reset!(sim.state.clock)
    return sim
end

"""
    $SIGNATURES

Advance the simulation forward by one timestep.
"""
timestep!(sim::Simulation) = timestep!(sim, get_dt(sim.time_stepping))
function timestep!(sim::Simulation, dt)
    reset_tendencies!(sim.state)
    timestep!(sim.state, sim.model, sim.time_stepping, dt)
    tick_time!(sim.state.clock, dt)
end

"""
    $SIGNATURES

Run the simulation by `steps` or a `period` with `dt` timestep size (in seconds or Dates.Period).
"""
function run!(sim::Simulation; 
        steps::Union{Int, Nothing} = nothing,
        period::Union{Period, Nothing} = nothing,
        dt = get_dt(get_time_stepping(sim.model)))

    dt = convert_dt(dt)
    steps = get_steps(steps, period, dt)

    for _ in 1:steps
        timestep!(sim, dt)
    end

    return sim 
end

current_time(sim::Simulation) = sim.state.clock.time

convert_dt(dt::Real) = dt
convert_dt(dt::Period) = Second(dt).value

get_steps(steps::Nothing, period::Period, dt::Real) = div(Second(period).value, dt)
get_steps(steps::Int, period::Nothing, dt::Real) = steps
get_steps(steps::Nothing, period::Nothing, dt::Real) = throw(ArgumentError("either `steps` or `period` must be specified"))
get_steps(steps::Int, period::Period, dt::Real) = throw(ArgumentError("both `steps` and `period` cannot be specified"))
