"""
    $TYPEDEF

Represents a "integrator" for a simulation of a given `model`. `ModelIntegrator` consists of a
`clock`, a `model`, and an initialized `StateVariables` data structure, as well as any relevant
`inputs` provided by a corresponding `InputProvider`. 
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
end

# Outer constructor so that `ModelIntegrator` can be constructed with a timestepper type
# so that it can correctly subtype Oceananigans.AbstractModel
# TODO: Is this actually needed?
# TODO: For now just take the explicit one?
function ModelIntegrator(
        clock::Clock,
        model::AbstractModel{NF},
        inputs::InputSources,
        state::AbstractStateVariables,
        initializers::NamedTuple,
    ) where {NF}
    grid = get_grid(model)
    timestepper = get_timestepper(model)
    return ModelIntegrator{
        NF, typeof(architecture(grid)), typeof(grid), typeof(timestepper),
        typeof(model), typeof(state), typeof(initializers), typeof(inputs),
    }(clock, model, inputs, state, initializers)
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

Oceananigans.Simulations.timestepper(integrator::ModelIntegrator) = get_timestepper(integrator.model)

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
    # set inputs based on updated clock/state
    initialize!(integrator.state, integrator.inputs)
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
    $TYPEDSIGNATURES

Advance the model forward by one timestep with optional timestep size `Δt`. If `finalize = true`,
`compute_auxiliary!` is called after the time step in order to update the values of auxiliary/diagnostic
variables.
"""
timestep!(integrator::ModelIntegrator; finalize = true) = timestep!(integrator, default_dt(get_timestepper(integrator.model)); finalize)
function timestep!(integrator::ModelIntegrator, Δt; finalize = true)
    timestep!(integrator, get_timestepper(integrator.model), convert_dt(Δt))
    if finalize
        compute_auxiliary!(integrator.state, integrator.model)
    end
    return nothing
end

"""
    timestep!(integrator::ModelIntegrator, timestepper::AbstractTimeStepper, Δt)

Advance the model forward by one timestep of size `Δt` using a single `timestepper`, which integrates
*all* prognostic variables. The variable names are forwarded to the corresponding
`timestep!(integrator, timestepper, Δt, names)` method and the clock is advanced once for the whole step.
"""
function timestep!(integrator::ModelIntegrator, timestepper::AbstractTimeStepper, Δt)
    # a single time stepper integrates all prognostic variables
    names = prognostic_names(integrator.state)
    isempty(names) || timestep!(integrator, timestepper, Δt, names)
    # advance the clock once for the entire step
    tick!(integrator.state.clock, Δt)
    return nothing
end

"""
    default_dt(integrator::ModelIntegrator)

Return the default timestep size for the given `integrator`, taken from its model's `timestepper`.
"""
default_dt(integrator::ModelIntegrator) = default_dt(get_timestepper(integrator.model))

"""
    $TYPEDSIGNATURES

Creates and initializes a `ModelIntegrator` for the given `model` with input variables populated by
the given `inputs`. This method allocates all necessary `Field`s for the state variables and subsequently calls
`initialize!(::ModelIntegrator)`.

The `timestepper_classes` keyword, a `NamedTuple` of `varname => class`, overrides the default timestepper class
(`:explicit`/`:implicit`) of the named prognostic variables (see [`prognostic`](@ref) for how defaults are declared).

Note that this method is **not type stable** and thus should not be called from Enzyme `autodiff`. To reinitialize
the model for an existing `state`, use `initialize!(state, model)`.
"""
function initialize(
        model::AbstractModel{NF},
        inputs::InputSource...;
        clock::Clock = Clock(time = zero(NF)),
        timestepper_classes = (;),
        boundary_conditions = (;),
        initializers = (;),
        fields = (;)
    ) where {NF}
    inputs = InputSources(inputs...)
    input_vars = variables(inputs)
    state = StateVariables(model; clock, timestepper_classes, boundary_conditions, fields, input_variables = input_vars)
    integrator = ModelIntegrator(clock, model, inputs, state, initializers)
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
    tsstr = get_timestepper(integrator.model)
    println(io, "Integrator of $modelstr with timestepper $tsstr")
    println(io, "├── Current time: $(current_time(integrator))")
    return println(io, "├── $statestr")
    # TODO: add more information?
end
