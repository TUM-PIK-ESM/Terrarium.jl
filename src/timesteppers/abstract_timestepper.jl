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
    $SIGNATURES

Allocate the time stepper cache for `timestepper` against the given `state`.
"""
initialize_timestepper_cache(timestepper::AbstractTimeStepper, state) = initialize(timestepper, state)

"""
    $SIGNATURES

Return the names of the prognostic variables in `state` that are integrated by the timestepper filling the
given `class` (`:explicit` or `:implicit`); i.e. those declared with `timestepper = class`. The names and
their classes are read from the `state` type parameters, which keeps this type stable.
"""
prognostic_names(state::StateVariables, class::Symbol) = prognostic_names(state, Val(class))

# Select, at compile time, the names of the prognostic variables whose timestepper class matches `class`.
# Operating on the `timestepper_classes` type parameters of `StateVariables` keeps this type stable.
@generated function prognostic_names(
        ::StateVariables{NF, prognames, closurenames, timestepper_classes}, ::Val{class}
    ) where {NF, prognames, closurenames, timestepper_classes, class}
    selected = Symbol[name for (name, c) in zip(prognames, timestepper_classes) if c === class]
    return Expr(:tuple, map(QuoteNode, selected)...)
end

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

Initialize and return the time stepping state cache for the given time stepper.
. Allocate and return a `NamedTuple` of intermediate fields/state required by the given
`timestepper`. Time steppers that do not require any cache can fall back to the default 
implementation, which returns an empty `NamedTuple`.
"""
initialize(timestepper::AbstractTimeStepper, state) = (;)

"""
    get_cache(state, timestepper::AbstractTimeStepper)

Return the cache associated with `timestepper` (retrieved from `state.timestepper_cache`). Falls back to
`nothing` for time steppers that do not define a cache.
"""
get_cache(state, ::AbstractTimeStepper) = nothing

"""
    explicit_cache(timestepper_cache)
    implicit_cache(timestepper_cache)

Select the cache belonging to the explicit/implicit (sub-)stepper from a `state.timestepper_cache`. For a
single timestepper the cache is stored bare and returned as-is; for an [`IMEX`](@ref) timestepper the
per-class sub-caches are selected by the IMEX-specific methods in `imex.jl`.
"""
explicit_cache(cache) = cache
implicit_cache(cache) = cache

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
