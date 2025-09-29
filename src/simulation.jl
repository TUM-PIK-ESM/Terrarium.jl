abstract type AbstractSimulation end

"""
    $TYPEDEF

Represents the state of a "simulation" for a given `model`. A simulation
is here defined as a particular choice of `model` in conjunction with
values for all state variables (including the `Clock`) and any necessary
state caches for the time-stepper.
"""
struct Simulation{
    Model<:AbstractModel,
    TimeStepperCache<:AbstractTimeStepperCache,
    StateVars<:AbstractStateVariables,
    Inputs<:InputProvider
} <: AbstractSimulation
    "The type of model used for the simulation."
    model::Model

    "Time stepper state cache."
    cache::TimeStepperCache

    "Input provider type"
    inputs::Inputs

    "Collection of all state variables defined on the simulation `model`."
    state::StateVars
end

"""
    $TYPEDEF

Creates and initializes a `Simulation` for the given `model` with the given `clock` state.
This method allocates all necessary `Field`s for the state variables and calls `initialize!(sim)`.
Note that this method is **not type stable** and should not be called in an Enzyme `autodiff` call.
"""
function initialize(model::AbstractModel{NF}, inputs::InputProvider; clock::Clock=Clock(time=zero(NF))) where {NF}
    state = StateVariables(model, clock, inputs.fields)
    time_stepping_cache = initialize(model, get_time_stepping(model))
    sim = Simulation(model, time_stepping_cache, inputs, state)
    initialize!(sim)
    return sim
end

# Convenience dispatch that constructs an InputProvider from zero or more input sources
function initialize(model::AbstractModel{NF}, inputs::InputSource...; clock::Clock=Clock(time=zero(NF))) where {NF}
    provider = InputProvider(get_grid(model), inputs...)
    return initialize(model, provider; clock)
end

"""
    $TYPEDEF

Resets the simulation `clock` and calls `initialize!(state, model)` on the underlying model which
should reset all state variables to their values as defiend by the model initializer.
"""
function initialize!(sim::Simulation)
    # TODO: reset other variables too?
    reset_tendencies!(sim.state)
    reset!(sim.state.clock)
    update_inputs!(sim.inputs, sim.state.clock)
    initialize!(sim.state, sim.model)
    return sim
end

"""
    $SIGNATURES

Advance the simulation forward by one timestep with optional timestep size `dt`.
"""
timestep!(sim::Simulation) = timestep!(sim, get_dt(get_time_stepping(sim.model)))
function timestep!(sim::Simulation, dt)
    reset_tendencies!(sim.state)
    update_inputs!(sim.inputs, sim.state.clock)
    timestep!(sim.state, sim.model, dt)
    tick_time!(sim.state.clock, dt)
end

"""
    $SIGNATURES

Run the simulation by `steps` or a `period` with `dt` timestep size (in seconds or Dates.Period).
"""
function run!(sim::Simulation;
    steps::Union{Int, Nothing} = nothing,
    period::Union{Period, Nothing} = nothing,
    dt = get_dt(get_time_stepping(sim.model))
)
    dt = convert_dt(dt)
    steps = get_steps(steps, period, dt)

    for _ in 1:steps
        timestep!(sim, dt)
    end

    return sim 
end

current_time(sim::Simulation) = sim.state.clock.time

get_steps(steps::Nothing, period::Period, dt::Real) = div(Second(period).value, dt)
get_steps(steps::Int, period::Nothing, dt::Real) = steps
get_steps(steps::Nothing, period::Nothing, dt::Real) = throw(ArgumentError("either `steps` or `period` must be specified"))
get_steps(steps::Int, period::Period, dt::Real) = throw(ArgumentError("both `steps` and `period` cannot be specified"))
