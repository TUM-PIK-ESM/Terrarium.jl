"""
Base type for time-stepper state caches.
"""
abstract type AbstractTimeStepperCache{NF} end

# AbstractTimeStepper

"""
Base type for time steppers.
"""
abstract type AbstractTimeStepper{NF} end

"""
    get_dt(timestepper::AbstractTimeStepper)

Get the current timestep size for the time stepper.
"""
function get_dt end

"""
    is_adaptive(timestepper::AbstractTimeStepper)

Returns `true` if the given time stepper is adaptive, false otherwise.
"""
function is_adaptive end

"""
    timestep!(state, model::AbstractModel, timestepper::AbstractTimeStepper, dt)

Advance the model state by one time step, or by `dt` units of time.
"""
function timestep! end

"""
    $SIGNATURES

Evaluates an explicit update `u += ∂u∂t*dt` for all prognostic variables `u` in `state`.
For prognostic variables with closures, `u` is the corresponding closure variable and `invclosure!`
is called after the update. If necessary, `explicit_step!` and/or `explicit_step_kernel!` can be
overridden for subtypes `AbstractTimeStepper` to implement specialized behavior.
"""
function explicit_step!(state, model::AbstractModel, timestepper::AbstractTimeStepper, dt, args...)
    # Get fields
    progfields = prognostic_fields(state)
    tendencies = tendency_fields(state)
    closurefields = closure_fields(state)
    # Evaluate step
    explicit_step!(progfields, closurefields, tendencies, get_grid(model), timestepper, dt, args...)
    # Update clock
    tick_time!(state.clock, dt)
    # Evaluate inverse closure relations (i.e. map from closure fields to prognostic fields)
    for name in keys(state.closures)
        invclosure!(state, model, state.closures[name])
    end
end

# Case 1: Named tuple (namespace) of prognostic/tendency variables; note that the names of the prognostic
# and tendnecy variables must match!
function explicit_step!(
    progfields::NamedTuple{names},
    closurefields::NamedTuple{cnames},
    tendencies::NamedTuple{names},
    grid::AbstractLandGrid,
    timestepper::AbstractTimeStepper,
    dt,
    args...
) where {names, cnames}
    for name in names
        if name ∈ cnames
            # if `name` is defined in `closurefields`, then it is either a prognostic field with
            # a closure relation or it is a namespace
            explicit_step!(progfields[name], closurefields[name], tendencies[name], grid, timestepper, dt, args...)
        else
            # if `name` is not in `closurefields`, then it must be a non-closure prognostic field
            explicit_step!(progfields[name], nothing, tendencies[name], grid, timestepper, dt, args...)
        end
    end
end

# Case 2a: non-closure prognostic field where tendency is applied directly
function explicit_step!(
    progfield::AbstractField{LX, LY, LZ},
    ::Nothing,
    tendency::AbstractField{LX, LY, LZ},
    grid::AbstractLandGrid,
    timestepper::AbstractTimeStepper,
    dt,
    args...
) where {LX, LY, LZ}
    launch!(
        grid,
        workspec(LX(), LY(), LZ()),
        explicit_step_kernel!,
        progfield,
        tendency,
        timestepper,
        dt,
        args...
    )
end

# Case 2b: prognostic field where the tendency is applied to the closure field
function explicit_step!(
    progfield::AbstractField{LX, LY, LZ},
    closurefield::AbstractField{LX, LY, LZ},
    tendency::AbstractField{LX, LY, LZ},
    grid::AbstractLandGrid,
    timestepper::AbstractTimeStepper,
    dt,
    args...
) where {LX, LY, LZ}
    launch!(
        grid,
        workspec(LX(), LY(), LZ()),
        explicit_step_kernel!,
        closurefield,
        tendency,
        timestepper,
        dt,
        args...,
    )
end

"""
    explicit_step_kernel!(field, tendency, timestepper::AbstractTimeStepper, dt)

Kernel that updates `field` based on the given `tendency` and `dt`. By default, this is
implemented as a simple Euler update `field += tendency*dt` which can serve as a building
block for more complex, multi-stage timesteppers. Where necessary, additional dispatches of
`explicit_step_kernel!` can be defined for specialized time-stepping schemes.
"""
@kernel function explicit_step_kernel!(field, tendency, ::AbstractTimeStepper, dt)
    i, j, k = @Index(Global, NTuple)

    @inbounds let
        u = field,
        ∂u∂t = tendency,
        dt = convert(eltype(∂u∂t), dt);
        u[i, j, k] += ∂u∂t[i, j, k] * dt
    end
end
