"""
    $TYPEDEF

Represents a "driver" for a simulation of a given `model`. `ModelDriver` consists of a
`clock`, a `model`, and an initialized `StateVariables` data structure, as well as a `stage`
for the timestepper and any relevant `inputs` provided by a corresponding `InputProvider`.
The `ModelDriver` implements the `Oceananigans.AbstractModel` interface and can thus be
treated as a "model" in `Oceananigans` `Simulation`s and output reading/writing utilities.
"""
struct ModelDriver{
    NF,
    Arch<:AbstractArchitecture,
    Grid<:AbstractLandGrid{NF, Arch},
    TimeStepper<:AbstractTimeStepper{NF},
    Model<:AbstractModel{NF, Grid},
    StateVars<:AbstractStateVariables,
    Inputs<:InputSources,
} <: Oceananigans.AbstractModel{TimeStepper, Arch}
    "The clock holding all information about the current timestep/iteration of a simulation"
    clock::Clock

    "Spatial discretization of the underlying `model`"
    grid::Grid

    "Underlying model evaluated by this driver"
    model::Model

    "Input sources"
    inputs::Inputs

    "Collection of all state variables defined on the simulation `model`"
    state::StateVars

    "Time stepper"
    timestepper::TimeStepper
end

# Oceananigans.AbstractModel interface

Base.time(driver::ModelDriver) = driver.clock.time

Base.eltype(::ModelDriver{NF}) where {NF} = NF

iteration(driver::ModelDriver) = driver.clock.iteration

architecture(driver::ModelDriver) = architecture(get_grid(driver.model))

timestepper(driver::ModelDriver) = driver.timestepper

update_state!(driver::ModelDriver; compute_tendencies = true) = update_state!(driver.state, driver.model, driver.inputs; compute_tendencies)

# for now, just forward Oceananigans.time_step! to timestep!
# consider renaming later...
time_step!(driver::ModelDriver, Δt; kwargs...) = timestep!(driver, Δt)

"""
    $TYPEDEF

Resets the simulation `clock` and calls `initialize!(state, model)` on the underlying model which
should reset all state variables to their values as defiend by the model initializer.
"""
function initialize!(driver::ModelDriver)
    # TODO: reset other variables too?
    reset_tendencies!(driver.state)
    reset!(driver.clock)
    update_inputs!(driver.state, driver.inputs)
    fill_halo_regions!(driver.state)
    initialize!(driver.state, driver.model)
    return driver
end

# Terrarium method interfaces

current_time(driver::ModelDriver) = driver.clock.time

get_fields(driver::ModelDriver, queries...) = get_fields(driver.state, queries...)

"""
    $SIGNATURES

Advance the model forward by one timestep with optional timestep size `Δt`. If `finalize = true`,
`compute_auxiliary!` is called after the time step in order to update the values of auxiliary/diagnostic
variables.
"""
timestep!(driver::ModelDriver; finalize = true) = timestep!(driver, default_dt(timestepper(driver)); finalize)
function timestep!(driver::ModelDriver, Δt; finalize = true)
    update_state!(driver, compute_tendencies = true)
    timestep!(driver.state, driver.timestepper, driver.model, driver.inputs, convert_dt(Δt))
    if finalize
        compute_auxiliary!(driver.state, driver.model)
    end
    return nothing
end

"""
    $SIGNATURES

Run the simulation for `steps` or a given time `period` with timestep size `Δt` (in seconds or Dates.Period).
"""
function run!(
    driver::ModelDriver;
    steps::Union{Int, Nothing} = nothing,
    period::Union{Period, Nothing} = nothing,
    Δt = default_dt(timestepper(driver))
)
    Δt = convert_dt(Δt)
    steps = get_steps(steps, period, Δt)

    for _ in 1:steps
        timestep!(driver, Δt, finalize=false)
    end

    # Update auxiliary variables for final timestep
    compute_auxiliary!(driver.state, driver.model)
    return driver 
end

"""
    $TYPEDEF

Creates and initializes a `ModelDriver` for the given `model` with the given `clock` state.
This method allocates all necessary `Field`s for the state variables and calls `initialize!(::ModelDriver)`.
Note that this method is **not type stable** and should not be called in an Enzyme `autodiff` call.
"""
function initialize(
    model::AbstractModel{NF},
    timestepper::AbstractTimeStepper;
    clock::Clock = Clock(time=zero(NF)),
    input_sources::InputSources = InputSources(),
    boundary_conditions = (;),
    fields = (;)
) where {NF}
    grid = get_grid(model)
    input_vars = variables(input_sources)
    state = initialize(model; clock, boundary_conditions, fields, external_variables=input_vars)
    initialized_timestepper = initialize(timestepper, model, state)
    driver = ModelDriver(clock, grid, model, input_sources, state, initialized_timestepper)
    initialize!(driver)
    return driver
end

get_steps(steps::Nothing, period::Period, Δt::Real) = div(Second(period).value, Δt)
get_steps(steps::Int, period::Nothing, Δt::Real) = steps
get_steps(steps::Nothing, period::Nothing, Δt::Real) = throw(ArgumentError("either `steps` or `period` must be specified"))
get_steps(steps::Int, period::Period, Δt::Real) = throw(ArgumentError("both `steps` and `period` cannot be specified"))

function Base.show(io::IO, driver::ModelDriver)
    modelstr = summary(driver.model)
    statestr = summary(driver.state)
    tsstr = summary(driver.timestepper)
    println(io, "Driver of $modelstr")
    println(io, "├── Current time: $(current_time(driver))")
    println(io, "├── Timestepper: $tsstr")
    println(io, "├── $statestr")
    # TODO: add more information?
end