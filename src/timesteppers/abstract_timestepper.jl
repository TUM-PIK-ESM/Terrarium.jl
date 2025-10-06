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
    default_dt(timestepper::AbstractTimeStepper)

Get the current timestep size for the time stepper.
"""
function default_dt end

"""
    is_adaptive(timestepper::AbstractTimeStepper)

Return `true` if the given time stepper is adaptive, false otherwise.
"""
function is_adaptive end

"""
    timestep!(state, model::AbstractModel, timestepper::AbstractTimeStepper, dt)

Advance the model state by one time step, or by `dt` units of time.
"""
function timestep! end

"""
    initialize(::AbstractTimeStepper, prognostic_fields, closure_fields, tendencies) where {NF}

Initialize and return the time stepping state cache for the given time stepper.
"""
initialize(::AbstractTimeStepper, prognostic_fields, closure_fields, tendencies) = nothing

"""
    $SIGNATURES

Evaluate an explicit update `u += ∂u∂t*dt` for all prognostic fields and their corresponding
tendencies. By default, this is implemented as a simple Euler update `u += dudt*dt` which can
serve as a building block for more complex, multi-stage timesteppers. Where necessary,
additional dispatches of `explicit_step_kernel!` can be defined for specialized
time-stepping schemes.
"""
function explicit_step!(state, grid::AbstractLandGrid, timestepper::AbstractTimeStepper, dt)
    fastiterate(keys(state.prognostic)) do name
        # update prognostic or closure state variable
        if haskey(state.closures, name)
            closure = state.closures[name]
            cname = varname(getvar(closure))
            explicit_step!(state.auxiliary[cname], state.tendencies[cname], grid, timestepper, dt)
        else
            explicit_step!(state.prognostic[name], state.tendencies[name], grid, timestepper, dt)
        end

    end
    fastiterate(state.namespaces) do ns
        explicit_step!(state, grid, timestepper, dt)
    end
    return nothing
end

function explicit_step!(
    progfield::AbstractField{LX, LY, LZ},
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

@kernel function explicit_step_kernel!(field, tendency, ::AbstractTimeStepper, dt)
    i, j, k = @index(Global, NTuple)
    u = field
    ∂u∂t = tendency
    @inbounds let dt = convert(eltype(tendency), dt);
        u[i, j, k] += ∂u∂t[i, j, k] * dt
    end
end
