abstract type AbstractSimulation end

struct Simulation{
    Model<:AbstractModel,
    StateVars<:AbstractStateVariables,
    ClockType<:Clock,
} <: AbstractSimulation
    model::Model
    
    state::StateVars

    clock::ClockType
end

function initialize!(sim::Simulation)
    # TODO: reset other variables too?
    reset_tendencies!(sim.state)
    initialize!(sim.state, sim.model)
    reset!(sim.clock)
end

timestep!(sim::Simulation) = timestep!(sim, get_dt(get_time_stepping(sim.model)))
function timestep!(sim::Simulation, dt)
    timestepper = get_time_stepping(sim.model)
    reset_tendencies!(sim.state)
    timestep!(sim.state, sim.model, timestepper, dt)
    tick_time!(sim.clock, dt)
end
