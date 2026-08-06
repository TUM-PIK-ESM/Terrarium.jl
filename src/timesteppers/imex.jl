"""
    $TYPEDEF

Base type for implicit-explicit (IMEX) time steppers. An `AbstractIMEX` integrates each prognostic variable
with one of two sub-steppers depending on its [`timestepping`](@ref) class: variables of class [`Explicit`](@ref)
are stepped by the explicit sub-stepper and those of class [`Implicit`](@ref) by the implicit sub-stepper.

Concrete subtypes (e.g. [`IMEX`](@ref)) must provide the explicit and implicit sub-steppers via
[`explicit_timestepper`](@ref) and [`implicit_timestepper`](@ref); most other behavior is defined here against
`AbstractIMEX`.
"""
abstract type AbstractIMEX{NF} <: AbstractTimeStepper{NF} end

"""
    explicit_timestepper(imex::AbstractIMEX)

Return the explicit sub-stepper of the [`AbstractIMEX`](@ref) timestepper `imex`.
"""
@inline explicit_timestepper(imex::AbstractIMEX) = imex.explicit

"""
    implicit_timestepper(imex::AbstractIMEX)

Return the implicit sub-stepper of the [`AbstractIMEX`](@ref) timestepper `imex`.
"""
@inline implicit_timestepper(imex::AbstractIMEX) = imex.implicit

"""
    $TYPEDEF

Implicit-explicit (IMEX) time stepper that integrates each prognostic variable with one of two sub-steppers
depending on its [`timestepping`](@ref) class: variables of class [`Explicit`](@ref) are stepped by `explicit`,
and those of class [`Implicit`](@ref) by `implicit`.

Each variable's class is resolved from `timestepping(var, model, imex)`, which defaults to `Explicit()` for
all variables. To integrate selected variables implicitly, specialize [`timestepping`](@ref) on the relevant
variable and/or model types together with the IMEX timestepper. The resolved per-variable classes are stored
in the [`IMEXCache`](@ref) type parameter and used to route variables at each step.

Properties:
$(TYPEDFIELDS)
"""
struct IMEX{
        NF,
        E <: AbstractTimeStepper{NF},
        I <: AbstractTimeStepper{NF},
    } <: AbstractIMEX{NF}
    "Sub-stepper for prognostic variables of class `Explicit` (should have `timestepping(explicit) == Explicit()`)"
    explicit::E

    "Sub-stepper for prognostic variables of class `Implicit` (should have `timestepping(implicit) == Implicit()`)"
    implicit::I
end

"""
    IMEX(explicit, implicit)
    IMEX(; explicit, implicit)

Construct an [`IMEX`](@ref) time stepper from an `explicit` and an `implicit` sub-stepper (which must share
the same numerical type `NF`). Which prognostic variables are integrated implicitly is controlled by
specializing [`timestepping`](@ref).
"""
IMEX(; explicit, implicit) = IMEX(explicit, implicit)

default_dt(imex::AbstractIMEX) = default_dt(explicit_timestepper(imex))

is_adaptive(imex::AbstractIMEX) = is_adaptive(explicit_timestepper(imex)) || is_adaptive(implicit_timestepper(imex))

"""
    $TYPEDEF

Cache for an [`AbstractIMEX`](@ref) time stepper. Holds each sub-stepper's own cache; the resolved per-variable
timestepping classes are stored as the leading type parameter `classes` (a tuple of [`Explicit`](@ref)/[`Implicit`](@ref)
instances, in prognostic-variable order) so that routing each variable to its sub-stepper is type stable.

Properties:
$(TYPEDFIELDS)
"""
struct IMEXCache{
        classes,
        NF,
        EC <: AbstractTimeStepperCache{NF},
        IC <: AbstractTimeStepperCache{NF},
    } <: AbstractTimeStepperCache{NF}
    "Cache for the explicit sub-stepper"
    explicit::EC

    "Cache for the implicit sub-stepper"
    implicit::IC
end

# Construct an IMEXCache with the resolved `classes` given as the leading type parameter.
function IMEXCache{classes}(explicit::EC, implicit::IC) where {classes, NF, EC <: AbstractTimeStepperCache{NF}, IC <: AbstractTimeStepperCache{NF}}
    return IMEXCache{classes, NF, EC, IC}(explicit, implicit)
end

function Adapt.adapt_structure(to, cache::IMEXCache{classes}) where {classes}
    return IMEXCache{classes}(Adapt.adapt_structure(to, cache.explicit), Adapt.adapt_structure(to, cache.implicit))
end

"""
    timestepping(cache::IMEXCache)

Return the resolved timestepping class of every prognostic variable (a tuple of [`Explicit`](@ref)/[`Implicit`](@ref)
instances, in prognostic-variable order) held in the [`IMEXCache`](@ref) type parameter.
"""
timestepping(::IMEXCache{classes}) where {classes} = classes

# Build the IMEX cache: allocate each sub-stepper's cache and resolve every prognostic variable's class
# (via `timestepping(var, model, imex)`), stored as the cache's type parameter.
function initialize(imex::AbstractIMEX, state, progvars, model)
    classes = resolve_timestepping(imex, progvars, model)
    return IMEXCache{classes}(initialize(explicit_timestepper(imex), state), initialize(implicit_timestepper(imex), state))
end

# A sub-stepper of an IMEX reads its own slice of the cache, selected by its `timestepping` trait.
get_cache(cache::IMEXCache, timestepper::AbstractTimeStepper) = get_cache(cache, timestepping(timestepper))
get_cache(cache::IMEXCache, ::Explicit) = cache.explicit
get_cache(cache::IMEXCache, ::Implicit) = cache.implicit

"""
    $SIGNATURES

Resolve the timestepping class of each prognostic variable in `progvars` (a `NamedTuple` of
[`PrognosticVariable`](@ref)s) for the IMEX timestepper `imex` and owning `model`, returned as a tuple of
[`Explicit`](@ref)/[`Implicit`](@ref) instances in the same order as `progvars`. Each variable's class is
given by [`timestepping`](@ref)`(var, model, imex)`.
"""
function resolve_timestepping(imex::AbstractIMEX, progvars::NamedTuple, model)
    return map(var -> timestepping(var, model, imex), values(progvars))
end

"""
    timestep!(integrator::ModelIntegrator, timestepper::AbstractIMEX, Δt)

Advance the model forward by one timestep of size `Δt` using an [`AbstractIMEX`](@ref) timestepper. Each
prognostic variable is routed to the explicit or implicit sub-stepper according to the resolved classes stored
in the [`IMEXCache`](@ref) type; each sub-stepper fetches its own sub-cache via [`get_cache`](@ref). The clock
is advanced once for the whole step.
"""
function timestep!(integrator::ModelIntegrator, timestepper::AbstractIMEX, Δt)
    cache = integrator.state.timestepper_cache
    # type-stable split of the prognostic variables by their resolved class
    explicit_names = prognostic_names(integrator.state, cache, Explicit())
    implicit_names = prognostic_names(integrator.state, cache, Implicit())
    isempty(explicit_names) || timestep!(integrator, explicit_timestepper(timestepper), Δt, explicit_names)
    isempty(implicit_names) || timestep!(integrator, implicit_timestepper(timestepper), Δt, implicit_names)
    # advance the clock once for the entire step
    tick!(integrator.state.clock, Δt)
    return nothing
end

# Select, at compile time, the prognostic variable names whose resolved class matches `class`, by zipping
# the state's prognostic names with the IMEX cache's resolved class list (both in prognostic-variable order).
@generated function prognostic_names(::StateVariables{NF, prognames}, ::IMEXCache{classes}, ::Class) where {NF, prognames, classes, Class}
    selected = Symbol[name for (name, c) in zip(prognames, classes) if c isa Class]
    return Expr(:tuple, map(QuoteNode, selected)...)
end
