abstract type AbstractSimulation end

struct Simulation{
    Model<:AbstractModel,
    StateVars<:AbstractStateVariables,
} <: AbstractSimulation
    "The type of model used for the simulation."
    model::Model

    "Collection of all state variables defined on the simulation `model`."
    state::StateVars
end

function initialize(model::AbstractModel; clock::Clock=Clock(time=0.0))
    state = StateVariables(model, clock)
    sim = Simulation(model, state)
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
timestep!(sim::Simulation) = timestep!(sim, get_dt(get_time_stepping(sim.model)))
function timestep!(sim::Simulation, dt)
    reset_tendencies!(sim.state)
    timestep!(sim.state, sim.model, dt)
    tick_time!(sim.state.clock, dt)
end
