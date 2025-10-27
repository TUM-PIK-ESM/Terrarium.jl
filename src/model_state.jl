"""
    $TYPEDEF

Represents the state of a "simulation" for a given `model`. `ModelState` consists of a
`clock`, a `model`, and an initialized `StateVariables` data structure, as well as a `cache`
for the timestepper and any relevant `inputs` provided by a corresponding `InputProvider`.
The `ModelState` implements the `Oceananigans.AbstractModel` interface and can thus be
treated as a "model" in `Oceananigans` `Simulation`s and output reading/writing utilities.
"""
struct ModelState{
    NF,
    Arch<:AbstractArchitecture,
    Grid<:AbstractLandGrid{NF, Arch},
    TimeStepper<:AbstractTimeStepper{NF},
    Model<:AbstractModel{NF, Grid},
    StateVars<:AbstractStateVariables,
    Inputs<:InputProvider
} <: Oceananigans.AbstractModel{TimeStepper, Arch}
    "The clock holding all information about the current timestep/iteration of a simulation"
    clock::Clock

    "Spatial discretization of the underlying `model`"
    grid::Grid

    "The type of model used for the simulation."
    model::Model

    "Input provider type"
    inputs::Inputs

    "Collection of all state variables defined on the simulation `model`."
    state::StateVars

    "Time stepper"
    timestepper::TimeStepper
end

# Oceananigans.AbstractModel interface

Base.time(state::ModelState) = state.clock.time

Base.eltype(::ModelState{NF}) where {NF} = NF

iteration(state::ModelState) = state.clock.iteration

architecture(state::ModelState) = architecture(get_grid(state.model))

get_time_stepping(state::ModelState) = state.timestepper

function update_state!(state::ModelState; compute_tendencies = true)
    reset_tendencies!(state.state)
    update_inputs!(state.inputs, state.clock)
    compute_auxiliary!(state.state, state.model)
    if compute_tendencies
        compute_tendencies!(state.state, state.model)
    end
end

# for now, just forward Oceananigans.time_step! to timestep!
# consider renaming later...
time_step!(state::ModelState, Δt; kwargs...) = timestep!(state, Δt)

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

get_fields(state::ModelState, queries...) = get_fields(state.state, queries...)

"""
    $SIGNATURES

Advance the model forward by one timestep with optional timestep size `Δt`.
"""
timestep!(state::ModelState; finalize=true) = timestep!(state, default_dt(get_time_stepping(state)); finalize)
function timestep!(state::ModelState, Δt; finalize=true)
    reset_tendencies!(state.state)
    update_inputs!(state.inputs, state.clock)
    timestep!(state.state, state.model, get_time_stepping(state), convert_dt(Δt))
    if finalize
        compute_auxiliary!(state.state, state.model)
    end
    return nothing
end

"""
    $SIGNATURES

Run the simulation by `steps` or a `period` with `Δt` timestep size (in seconds or Dates.Period).
"""
function run!(
    state::ModelState;
    steps::Union{Int, Nothing} = nothing,
    period::Union{Period, Nothing} = nothing,
    Δt = default_dt(get_time_stepping(state))
)
    Δt = convert_dt(Δt)
    steps = get_steps(steps, period, Δt)

    for _ in 1:steps
        timestep!(state, Δt, finalize=false)
    end

    # Update auxiliary variables for final timestep
    compute_auxiliary!(state.state, state.model)

    return state 
end

"""
    $TYPEDEF

Creates and initializes a `ModelState` for the given `model` with the given `clock` state.
This method allocates all necessary `Field`s for the state variables and calls `initialize!(::ModelState)`.
Note that this method is **not type stable** and should not be called in an Enzyme `autodiff` call.
"""
function initialize(model::AbstractModel{NF}, inputs::InputProvider; clock::Clock=Clock(time=zero(NF)), timestepper::Type=ForwardEuler) where {NF}
    statevars = StateVariables(model, clock, inputs.fields)
    state = ModelState(clock, get_grid(model), model, inputs, statevars, timestepper(statevars))
    initialize!(state)
    return state
end

# Convenience dispatch that constructs an InputProvider from zero or more input sources
function initialize(model::AbstractModel{NF}, inputs::InputSource...; clock::Clock=Clock(time=zero(NF))) where {NF}
    provider = InputProvider(get_grid(model), inputs...)
    return initialize(model, provider; clock)
end

get_steps(steps::Nothing, period::Period, Δt::Real) = div(Second(period).value, Δt)
get_steps(steps::Int, period::Nothing, Δt::Real) = steps
get_steps(steps::Nothing, period::Nothing, Δt::Real) = throw(ArgumentError("either `steps` or `period` must be specified"))
get_steps(steps::Int, period::Period, Δt::Real) = throw(ArgumentError("both `steps` and `period` cannot be specified"))

function Base.show(io::IO, mime::MIME"text/plain", state::ModelState)
    println(io, "ModelState of $(typeof(state.model))")
    println(io, "  current time: $(current_time(state))")
    # TODO: add more information?
end
