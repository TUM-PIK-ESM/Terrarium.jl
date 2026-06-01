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
    $SIGNATURES

Return the class (`:explicit` or `:implicit`) that the given `timestepper` fills within a model's
`timesteppers` `NamedTuple`. This is determined by the abstract supertype of the time stepper.
"""
timestepper_class(::AbstractExplicitTimestepper) = :explicit
timestepper_class(::AbstractImplicitTimestepper) = :implicit

"""
    default_dt(timestepper::AbstractTimeStepper)
    default_dt(timesteppers::NamedTuple)

Get the current timestep size for the time stepper. For a `NamedTuple` of timesteppers, the timestep
size of the `explicit` timestepper is returned.
"""
function default_dt end

default_dt(timesteppers::NamedTuple) = default_dt(timesteppers.explicit)

"""
    is_adaptive(timestepper::AbstractTimeStepper)
    is_adaptive(timesteppers::NamedTuple)

Return `true` if the given time stepper is adaptive, false otherwise. For a `NamedTuple` of timesteppers,
returns `true` if *any* of the timesteppers are adaptive.
"""
function is_adaptive end

is_adaptive(timesteppers::NamedTuple) = any(is_adaptive, values(timesteppers))

"""
    $SIGNATURES

Default `timesteppers` for models, consisting of a single `explicit` [`ForwardEuler`](@ref) time stepper.
"""
default_timesteppers(::Type{NF}) where {NF} = (; explicit = ForwardEuler(NF))

# TODO: currently at most one per class of timestepper, we might relax that in the future
"""
    to_timesteppers(::Type{NF}, timesteppers)

Normalize a single timestepper, a tuple/vector of timesteppers, or a `NamedTuple` into the canonical
`timesteppers` `NamedTuple` with (at most) the entries `explicit` and `implicit`. If no explicit timestepper 
is provided, the `explicit` entry defaults to [`ForwardEuler`](@ref).
"""
to_timesteppers(::Type{NF}, timestepper::AbstractTimeStepper) where {NF} = to_timesteppers(NF, (timestepper,))
function to_timesteppers(::Type{NF}, timesteppers::Union{Tuple, AbstractVector}) where {NF}
    explicit_steppers = filter(ts -> isa(ts, AbstractExplicitTimestepper), timesteppers)
    implicit_steppers = filter(ts -> isa(ts, AbstractImplicitTimestepper), timesteppers)
    length(explicit_steppers) <= 1 || throw(ArgumentError("at most one explicit timestepper can be specified, got $(length(explicit_steppers))"))
    length(implicit_steppers) <= 1 || throw(ArgumentError("at most one implicit timestepper can be specified, got $(length(implicit_steppers))"))
    explicit = isempty(explicit_steppers) ? ForwardEuler(NF) : first(explicit_steppers)
    return isempty(implicit_steppers) ? (; explicit) : (; explicit, implicit = first(implicit_steppers))
end
function to_timesteppers(::Type{NF}, timesteppers::NamedTuple) where {NF}
    all(in((:explicit, :implicit)), keys(timesteppers)) || throw(ArgumentError("`timesteppers` keys must be a subset of (:explicit, :implicit), got $(keys(timesteppers))"))
    explicit = get(timesteppers, :explicit, ForwardEuler(NF))
    return haskey(timesteppers, :implicit) ? (; explicit, implicit = timesteppers.implicit) : (; explicit)
end

"""
    get_timesteppers(model::AbstractModel)::NamedTuple

Return the `timesteppers` associated with the given `model`. All `AbstractModel`s are required to
define a `timesteppers` field.
"""
@inline get_timesteppers(model::AbstractModel) = model.timesteppers

"""
    $SIGNATURES

Return the names of the prognostic variables defined by `model` that are integrated by the timestepper
filling the given `class` (`:explicit` or `:implicit`); i.e. those declared with `timestepper = class`.
"""
function prognostic_names_for(model::AbstractModel, class::Symbol)
    progvars = prognostic_variables(model)
    selected = filter(var -> timestepper(var) === class, progvars)
    return map(varname, selected)
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

Return the cache associated with `timestepper` (stored in `state.timestepper_cache` under the
timestepper's class). Falls back to `nothing` for time steppers that do not define a cache.
"""
get_cache(state, ::AbstractTimeStepper) = nothing

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
