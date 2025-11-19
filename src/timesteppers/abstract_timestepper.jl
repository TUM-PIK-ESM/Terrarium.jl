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
    is_initialized(timestepper::AbstractTimeStepper)

Return `true` if the timestepper is initialized, `false` otherwise.
"""
function is_initialized end

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
    timestep!(state, timestepper::AbstractTimeStepper, model::AbstractModel, inputs::InputSources, Δt)

Advance prognostic variables by one time step based on the current state, or by `Δt` units of time.
"""
function timestep! end

"""
    initialize(::AbstractTimeStepper, model, state) where {NF}

Initialize and return the time stepping state cache for the given time stepper.
"""
initialize(timestepper::AbstractTimeStepper, model, state) = timestepper

"""
    $SIGNATURES

Evaluate an explicit update `u += ∂u∂t*Δt` for all prognostic fields and their corresponding
tendencies. By default, this is implemented as a simple Euler update `u += dudt*Δt` which can
serve as a building block for more complex, multi-stage timesteppers. Where necessary,
additional dispatches of `explicit_step_kernel!(field, tendency, ::AbstractLandGrid, ::TimeStepper, Δt)`
can be defined to implement more specialized time-stepping schemes.
"""
function explicit_step!(state, grid::AbstractLandGrid, timestepper::AbstractTimeStepper, Δt)
    @assert is_initialized(timestepper)
    fastiterate(keys(state.prognostic)) do name
        # update prognostic state variable
        explicit_step!(state.prognostic[name], state.tendencies[name], grid, timestepper, Δt)
    end
    fastiterate(state.namespaces) do ns
        explicit_step!(ns, grid, timestepper, Δt)
    end
    return nothing
end

"""
Accumulate `tendency*Δt` in the given prognostic `field`. This method can be overridden by specialized
timestepping schemes as needed.
"""
function explicit_step!(
        field::AbstractField{LX, LY, LZ},
        tendency::AbstractField{LX, LY, LZ},
        grid::AbstractLandGrid,
        timestepper::AbstractTimeStepper,
        Δt,
        args...
    ) where {LX, LY, LZ}
    return launch!(
        grid,
        :xyz,
        explicit_step_kernel_3D!,
        field,
        tendency,
        timestepper,
        Δt,
        args...
    )
end

function explicit_step!(
        field::AbstractField{LX, LY, Nothing},
        tendency::AbstractField{LX, LY, Nothing},
        grid::AbstractLandGrid,
        timestepper::AbstractTimeStepper,
        Δt,
        args...
    ) where {LX, LY}
    return launch!(
        grid,
        :xy,
        explicit_step_kernel_2D!,
        field,
        tendency,
        timestepper,
        Δt,
        args...
    )
end

@kernel function explicit_step_kernel_3D!(
        field,
        tendency,
        ::AbstractTimeStepper,
        Δt
    )
    i, j, k = @index(Global, NTuple)
    u = field
    ∂u∂t = tendency
    @inbounds let Δt = convert(eltype(tendency), Δt)
        u[i, j, k] += ∂u∂t[i, j, k] * Δt
    end
end

@kernel function explicit_step_kernel_2D!(
        field,
        tendency,
        ::AbstractTimeStepper,
        Δt
    )
    i, j = @index(Global, NTuple)
    u = field
    ∂u∂t = tendency
    @inbounds let Δt = convert(eltype(tendency), Δt)
        u[i, j, 1] += ∂u∂t[i, j] * Δt
    end
end
