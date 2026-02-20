"""
    $TYPEDEF

Represents a "integrator" for a simulation of a given `model`. `ModelIntegrator` consists of a
`clock`, a `model`, and an initialized `StateVariables` data structure, as well as a `stage`
for the timestepper and any relevant `inputs` provided by a corresponding `InputProvider`.
The `ModelIntegrator` implements the `Oceananigans.AbstractModel` interface and can thus be
treated as a "model" in `Oceananigans` `Simulation`s and output reading/writing utilities.
"""
struct ModelIntegrator{
        NF,
        Arch <: AbstractArchitecture,
        Grid <: AbstractLandGrid{NF, Arch},
        TimeStepper <: AbstractTimeStepper{NF},
        Model <: AbstractModel{NF, Grid},
        StateVars <: AbstractStateVariables,
        Inits <: NamedTuple,
        Inputs <: InputSources,
    } <: Oceananigans.AbstractModel{TimeStepper, Arch}
    "The clock holding all information about the current timestep/iteration of a simulation"
    clock::Clock

    "Underlying model evaluated by this integrator"
    model::Model

    "Input sources"
    inputs::Inputs

    "Collection of all state variables defined on the simulation `model`"
    state::StateVars

    "Optional named tuple of user-specified field initializers"
    initializers::Inits

    "Time stepper"
    timestepper::TimeStepper
end

# Oceananigans model interface

Base.time(integrator::ModelIntegrator) = integrator.clock.time

Base.eltype(::ModelIntegrator{NF}) where {NF} = NF

function Base.getproperty(integrator::ModelIntegrator, name::Symbol)
    # Temporary hack to make Oceananigans output writers play nicely with ModelIntegrator
    # TODO: Raise an issue in Oceananigans for better long-term solution
    if name == :grid
        model = getfield(integrator, :model)
        return get_field_grid(get_grid(model))
    else
        return getfield(integrator, name)
    end
end

Oceananigans.Solvers.iteration(integrator::ModelIntegrator) = integrator.clock.iteration

Oceananigans.Architectures.architecture(integrator::ModelIntegrator) = architecture(get_grid(integrator.model))

Oceananigans.TimeSteppers.update_state!(integrator::ModelIntegrator; compute_tendencies = true) = update_state!(integrator.state, integrator.model, integrator.inputs; compute_tendencies)

# for now, just forward Oceananigans.time_step! to timestep!
# consider renaming later...
Oceananigans.TimeSteppers.time_step!(integrator::ModelIntegrator, Δt; kwargs...) = timestep!(integrator, Δt)

Oceananigans.Simulations.timestepper(integrator::ModelIntegrator) = integrator.timestepper
"""
    $SIGNATURES

Run the simulation for `steps` or a given time `period` with timestep size `Δt` (in seconds or Dates.Period).
"""
function Oceananigans.Simulations.run!(
        integrator::ModelIntegrator;
        steps::Union{Int, Nothing} = nothing,
        period::Union{Period, Nothing} = nothing,
        Δt = default_dt(timestepper(integrator))
    )
    Δt = convert_dt(Δt)
    steps = get_steps(steps, period, Δt)

    for _ in 1:steps
        timestep!(integrator, Δt, finalize = false)
    end

    # Update auxiliary variables for final timestep
    compute_auxiliary!(integrator.state, integrator.model)
    return integrator
end

"""
    $TYPEDEF

Resets the simulation `clock` and calls `initialize!(state, model)` on the underlying model which
should reset all state variables to their values as defiend by the model initializer.
"""
function initialize!(integrator::ModelIntegrator)
    # reset state variables and clock
    reset!(integrator.state)
    reset!(integrator.clock)
    # set inptus based on updated clock/state
    update_inputs!(integrator.state, integrator.inputs)
    # fill halo regions
    fill_halo_regions!(integrator.state)
    # evaluate user-specified field initializers
    initialize!(integrator.state, integrator.initializers)
    # evaluate model initializer
    initialize!(integrator.state, integrator.model)
    return integrator
end

# Terrarium method interfaces

current_time(integrator::ModelIntegrator) = integrator.clock.time

get_fields(integrator::ModelIntegrator, queries...) = get_fields(integrator.state, queries...)

"""
    $SIGNATURES

Advance the model forward by one timestep with optional timestep size `Δt`. If `finalize = true`,
`compute_auxiliary!` is called after the time step in order to update the values of auxiliary/diagnostic
variables.
"""
timestep!(integrator::ModelIntegrator; finalize = true) = timestep!(integrator, default_dt(timestepper(integrator)); finalize)
function timestep!(integrator::ModelIntegrator, Δt; finalize = true)
    update_state!(integrator, compute_tendencies = true)
    timestep!(integrator.state, integrator.timestepper, integrator.model, integrator.inputs, convert_dt(Δt))
    if finalize
        compute_auxiliary!(integrator.state, integrator.model)
    end
    return nothing
end

"""
    $TYPEDEF

Creates and initializes a `ModelIntegrator` for the given `model` with the given `clock` state.
This method allocates all necessary `Field`s for the state variables and calls `initialize!(::ModelIntegrator)`.
Note that this method is **not type stable** and should not be called in an Enzyme `autodiff` call.
"""
function initialize(
        model::AbstractModel{NF},
        timestepper::AbstractTimeStepper,
        inputs::InputSource...;
        clock::Clock = Clock(time = zero(NF)),
        boundary_conditions = (;),
        initializers = (;),
        fields = (;)
    ) where {NF}
    inputs = InputSources(inputs...)
    input_vars = variables(inputs)
    state = initialize(model; clock, boundary_conditions, fields, input_variables = input_vars)
    initialized_timestepper = initialize(timestepper, model, state)
    integrator = ModelIntegrator(clock, model, inputs, state, initializers, initialized_timestepper)
    initialize!(integrator)
    return integrator
end

get_steps(steps::Nothing, period::Period, Δt::Real) = div(Second(period).value, Δt)
get_steps(steps::Int, period::Nothing, Δt::Real) = steps
get_steps(steps::Nothing, period::Nothing, Δt::Real) = throw(ArgumentError("either `steps` or `period` must be specified"))
get_steps(steps::Int, period::Period, Δt::Real) = throw(ArgumentError("both `steps` and `period` cannot be specified"))

function Base.show(io::IO, integrator::ModelIntegrator)
    modelstr = summary(integrator.model)
    statestr = summary(integrator.state)
    tsstr = summary(integrator.timestepper)
    println(io, "Integrator of $modelstr with $tsstr")
    println(io, "├── Current time: $(current_time(integrator))")
    return println(io, "├── $statestr")
    # TODO: add more information?
end
