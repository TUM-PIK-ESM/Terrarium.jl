"""
    $TYPEDEF

Represents the state of a "simulation" for a given `model`. `ModelState` consists of a
`clock`, a `model`, and an initialized `StateVariables` data structure, as well as a `cache`
for the timestepper and any relevant `inputs` provided by a corresponding `InputProvider`.
The `ModelState` implements the `Oceananigans.AbstractModel` interface and can thus be
provided as a "model" to all `Oceananigans` types.
"""
struct ModelState{
    NF,
    Arch<:AbstractArchitecture,
    Grid<:AbstractLandGrid{NF, Arch},
    TimeStepper<:AbstractTimeStepper,
    Model<:AbstractModel{NF, Grid, TimeStepper},
    TimeStepperCache<:Union{Nothing, AbstractTimeStepperCache},
    StateVars<:AbstractStateVariables,
    Inputs<:InputProvider
} <: Oceananigans.AbstractModel{TimeStepper, Arch}
    "The clock holding all information about the current timestep/iteration of a simulation"
    clock::Clock

    "The type of model used for the simulation."
    model::Model

    "Time stepper state cache."
    cache::TimeStepperCache

    "Input provider type"
    inputs::Inputs

    "Collection of all state variables defined on the simulation `model`."
    state::StateVars
end

# Oceananigans.AbstractModel interface

Base.time(state::ModelState) = state.clock.time

iteration(state::ModelState) = state.clock.iteration

architecture(state::ModelState) = architecture(get_grid(state.model))

timestepper(state::ModelState) = get_time_stepping(state.model)

# for now, just forward Oceananigans.time_step! to timestep!
# consider renaming later...
time_step!(state::ModelState, args...) = timestep!(state, args...)

"""
    $TYPEDEF

Resets the simulation `clock` and calls `initialize!(state, model)` on the underlying model which
should reset all state variables to their values as defiend by the model initializer.
"""
function initialize!(state::ModelState)
    # TODO: reset other variables too?
    reset_tendencies!(state.state)
    reset!(state.clock)
    update_inputs!(state.inputs, state.clock)
    initialize!(state.state, state.model)
    return state
end

# Terrarium method interfaces

current_time(state::ModelState) = state.clock.time

"""
    $SIGNATURES

Advance the model forward by one timestep with optional timestep size `dt`.
"""
timestep!(state::ModelState) = timestep!(state, get_dt(get_time_stepping(state.model)))
function timestep!(state::ModelState, dt)
    reset_tendencies!(state.state)
    update_inputs!(state.inputs, state.clock)
    timestep!(state.state, state.model, dt)
    tick_time!(state.clock, dt)
end

"""
    $SIGNATURES

Run the simulation by `steps` or a `period` with `dt` timestep size (in seconds or Dates.Period).
"""
function run!(
    state::ModelState;
    steps::Union{Int, Nothing} = nothing,
    period::Union{Period, Nothing} = nothing,
    dt = get_dt(get_time_stepping(state.model))
)
    dt = convert_dt(dt)
    steps = get_steps(steps, period, dt)

    for _ in 1:steps
        timestep!(state, dt)
    end

    return state 
end

"""
    $TYPEDEF

Creates and initializes a `ModelState` for the given `model` with the given `clock` state.
This method allocates all necessary `Field`s for the state variables and calls `initialize!(::ModelState)`.
Note that this method is **not type stable** and should not be called in an Enzyme `autodiff` call.
"""
function initialize(model::AbstractModel{NF}, inputs::InputProvider; clock::Clock=Clock(time=zero(NF))) where {NF}
    statevars = StateVariables(model, clock, inputs.fields)
    time_stepping_cache = initialize(model, get_time_stepping(model))
    state = ModelState(clock, model, time_stepping_cache, inputs, statevars)
    initialize!(state)
    return state
end

# Convenience dispatch that constructs an InputProvider from zero or more input sources
function initialize(model::AbstractModel{NF}, inputs::InputSource...; clock::Clock=Clock(time=zero(NF))) where {NF}
    provider = InputProvider(get_grid(model), inputs...)
    return initialize(model, provider; clock)
end

get_steps(steps::Nothing, period::Period, dt::Real) = div(Second(period).value, dt)
get_steps(steps::Int, period::Nothing, dt::Real) = steps
get_steps(steps::Nothing, period::Nothing, dt::Real) = throw(ArgumentError("either `steps` or `period` must be specified"))
get_steps(steps::Int, period::Period, dt::Real) = throw(ArgumentError("both `steps` and `period` cannot be specified"))

function Base.show(io::IO, mime::MIME"text/plain", state::ModelState)
    println(io, "$(typeof(state.model)) state")
    println(io, "  current time: $(current_time(state))")
    # TODO: add more information?
end
