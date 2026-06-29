"""
    $TYPEDEF

Base type for time stepper caches. Each [`AbstractTimeStepper`](@ref) allocates a corresponding
`AbstractTimeStepperCache` subtype (via [`initialize`](@ref)) to hold any working state it needs between
stages/steps.
"""
abstract type AbstractTimeStepperCache{NF} end

"""
    $TYPEDEF

Trivial cache for time steppers that require no working state (e.g. [`ForwardEuler`](@ref)).
"""
struct EmptyCache{NF} <: AbstractTimeStepperCache{NF} end

# AbstractTimeStepper

"""
Base type for time steppers.
"""
abstract type AbstractTimeStepper{NF} end

"""
    $TYPEDEF

Base type for *explicit* time steppers.
"""
abstract type AbstractExplicitTimestepper{NF} <: AbstractTimeStepper{NF} end

"""
    $TYPEDEF

Base type for *implicit* time steppers.
"""
abstract type AbstractImplicitTimestepper{NF} <: AbstractTimeStepper{NF} end

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
    $SIGNATURES

Default `timestepper` for models: a single `explicit` [`ForwardEuler`](@ref) time stepper.
"""
default_timestepper(::Type{NF}) where {NF} = ForwardEuler(NF)

"""
    get_timestepper(model::AbstractModel)::AbstractTimeStepper

Return the `timestepper` associated with the given `model`. All `AbstractModel`s are required to
define a `timestepper` field holding an [`AbstractTimeStepper`](@ref) (e.g. [`ForwardEuler`](@ref),
[`Heun`](@ref), or [`IMEX`](@ref)).
"""
@inline get_timestepper(model::AbstractModel) = model.timestepper

"""
    initialize(timestepper::AbstractTimeStepper, state, progvars)

Allocate the time stepper cache for `timestepper` against the given `state`. `progvars` is the named tuple
of prognostic variable metadata (needed e.g. by [`IMEX`](@ref))
"""
initialize(timestepper::AbstractTimeStepper, state, progvars) = initialize(timestepper, state)

"""
    timestep!(integrator::ModelIntegrator, timestepper::AbstractTimeStepper, Δt)

Advance prognostic variables of the `integrator` model by one time step based on the current state, or by `Δt` units of time.
"""
function timestep! end

"""
    timestep!(state, model::AbstractModel, timestepper::AbstractTimeStepper, Δt)

Apply any necessary corrections or model-specific time stepping logic after applying `timestepper` to the prognostic state
variables defined by `model`.
"""
timestep!(state, model::AbstractModel, timestepper::AbstractTimeStepper, Δt) = nothing

"""
    initialize(::AbstractTimeStepper, state)

Initialize and return the [`AbstractTimeStepperCache`](@ref) holding any intermediate fields/state required by the
given `timestepper`. Time steppers that need no working state fall back to the default implementation,
which returns an [`EmptyCache`](@ref).
"""
initialize(timestepper::AbstractTimeStepper{NF}, state) where {NF} = EmptyCache{NF}()

"""
    get_cache(cache::AbstractTimeStepperCache, timestepper::AbstractTimeStepper)

Return the working cache for `timestepper` given the model's `state.timestepper_cache`. For a single
timestepper this is the cache itself; an [`IMEX`](@ref) cache returns the sub-cache matching the
timestepper's class.
"""
get_cache(cache::AbstractTimeStepperCache, ::AbstractTimeStepper) = cache

"""
    $SIGNATURES

Evaluate an explicit update `u += ∂u∂t*Δt` for the prognostic fields of `state` listed in `names` and
their corresponding tendencies. By default, this is implemented as a simple Euler update `u += dudt*Δt`
which can serve as a building block for more complex, multi-stage timesteppers. Where necessary,
additional dispatches of `explicit_step_kernel!(field, tendency, ::AbstractLandGrid, ::TimeStepper, Δt)`
can be defined to implement more specialized time-stepping schemes.
"""
function explicit_step!(state, grid::AbstractLandGrid, timestepper::AbstractTimeStepper, Δt, names::Tuple{Vararg{Symbol}})
    fastiterate(names) do name
        # apply flux BCs, if present
        compute_z_bcs!(state.tendencies[name], state.prognostic[name], grid, state)
        # debug site post-BC
        debugsite!(explicit_step!, state.tendencies[name], name)
        # update prognostic state variable
        explicit_step!(state.prognostic[name], state.tendencies[name], grid, timestepper, Δt)
        # debug site post-step
        debugsite!(explicit_step!, state.prognostic[name], name)
    end
    fastiterate(state.namespaces) do ns
        explicit_step!(ns, grid, timestepper, Δt, prognostic_names(ns))
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
    launch!(
        grid, XYZ, explicit_step_xyz_kernel!,
        field, tendency, timestepper, Δt, args...
    )
    return nothing
end

function explicit_step!(
        field::AbstractField{LX, LY, Nothing},
        tendency::AbstractField{LX, LY, Nothing},
        grid::AbstractLandGrid,
        timestepper::AbstractTimeStepper,
        Δt,
        args...
    ) where {LX, LY}
    launch!(
        grid, XY, explicit_step_xy_kernel!,
        field, tendency, timestepper, Δt, args...
    )
    return nothing
end

@kernel function explicit_step_xyz_kernel!(
        field,
        grid,
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

@kernel function explicit_step_xy_kernel!(
        field,
        grid,
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

# Default debug hooks
@inline debughook!(::typeof(explicit_step!), field, name) = checkfinite!(field, name)
