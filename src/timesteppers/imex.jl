"""
    $TYPEDEF

Implicit-explicit (IMEX) time stepper that integrates each prognostic variable with one of two
sub-steppers depending on its timestepper class: variables of class `:explicit` are stepped by `explicit`,
and those of class `:implicit` are stepped by `implicit`.

Each variable's class defaults to the one declared on its [`PrognosticVariable`](@ref) (see
[`prognostic`](@ref)) and can be overridden per-variable via the `timestepper_classes` field, a `NamedTuple`
of `varname => class`. The resolved per-variable classes are stored in the [`IMEXCache`](@ref) type
parameter and used to route variables at each step.

Properties:
$(TYPEDFIELDS)
"""
struct IMEX{
        NF,
        E <: AbstractExplicitTimestepper{NF},
        I <: AbstractImplicitTimestepper{NF},
        TC <: NamedTuple,
    } <: AbstractTimeStepper{NF}
    "Sub-stepper for prognostic variables of class `:explicit`"
    explicit::E

    "Sub-stepper for prognostic variables of class `:implicit`"
    implicit::I

    "Per-variable timestepper class overrides (`varname => :explicit`/`:implicit`)"
    timestepper_classes::TC
end

"""
    IMEX(explicit, implicit; timestepper_classes = (;))
    IMEX(; explicit, implicit, timestepper_classes = (;))

Construct an [`IMEX`](@ref) time stepper from an `explicit` and an `implicit` sub-stepper (which must share
the same numerical type). `timestepper_classes` is a `NamedTuple` of `varname => :explicit`/`:implicit` that
overrides the default class (declared on each [`PrognosticVariable`](@ref)) for the named variables.
"""
function IMEX(explicit::AbstractExplicitTimestepper{NF}, implicit::AbstractImplicitTimestepper{NF}; timestepper_classes::NamedTuple = (;)) where {NF}
    return IMEX{NF, typeof(explicit), typeof(implicit), typeof(timestepper_classes)}(explicit, implicit, timestepper_classes)
end
IMEX(; explicit, implicit, timestepper_classes = (;)) = IMEX(explicit, implicit; timestepper_classes)

default_dt(imex::IMEX) = default_dt(imex.explicit)

is_adaptive(imex::IMEX) = is_adaptive(imex.explicit) || is_adaptive(imex.implicit)

"""
    $TYPEDEF

Cache for an [`IMEX`](@ref) time stepper. Holds each sub-stepper's own cache; the resolved per-variable
timestepper classes are stored as the leading type parameter `classes` (a tuple of `:explicit`/`:implicit`
symbols, in prognostic-variable order) so that routing each variable to its sub-stepper is type stable.

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

"""
    timestepper_classes(cache::IMEXCache)

Return the resolved class of every prognostic variable (a tuple of `:explicit`/`:implicit` symbols, in
prognostic-variable order) held in the [`IMEXCache`](@ref) type parameter.
"""
timestepper_classes(::IMEXCache{classes}) where {classes} = classes

# Build the IMEX cache: allocate each sub-stepper's cache and resolve every prognostic variable's class
# (default from `progvars`, overridden by the IMEX overrides), stored as the cache's type parameter.
function initialize(imex::IMEX, state, progvars)
    classes = resolve_timestepper_classes(progvars, imex.timestepper_classes)
    return IMEXCache{classes}(initialize(imex.explicit, state), initialize(imex.implicit, state))
end

# A sub-stepper of an IMEX reads its own slice of the cache, selected by its class.
get_cache(cache::IMEXCache, ::AbstractExplicitTimestepper) = cache.explicit
get_cache(cache::IMEXCache, ::AbstractImplicitTimestepper) = cache.implicit

"""
    $SIGNATURES

Resolve the timestepper class of each prognostic variable in `progvars` (a `NamedTuple` of
[`PrognosticVariable`](@ref)s), returned as a tuple of `:explicit`/`:implicit` symbols in the same order as
`progvars`. Each variable's default class is the one it was declared with (see [`prognostic`](@ref)); the
`overrides` `NamedTuple` replaces the class for exactly the named variables. Unknown override keys error.
"""
function resolve_timestepper_classes(progvars::NamedTuple, overrides::NamedTuple)
    for key in keys(overrides)
        key in keys(progvars) || throw(ArgumentError(
            "`timestepper_classes` has unknown key :$key; expected a prognostic variable in $(keys(progvars))"
        ))
    end
    return map(values(progvars)) do var
        class = get(overrides, varname(var), timestepper(var))
        class in (:explicit, :implicit) || throw(ArgumentError(
            "timestepper class for prognostic variable :$(varname(var)) must be :explicit or :implicit, got :$class"
        ))
        return class
    end
end

"""
    timestep!(integrator::ModelIntegrator, timestepper::IMEX, Δt)

Advance the model forward by one timestep of size `Δt` using an [`IMEX`](@ref) timestepper. Each prognostic
variable is routed to the explicit or implicit sub-stepper according to the resolved classes stored in the
[`IMEXCache`](@ref) type; each sub-stepper fetches its own sub-cache via [`get_cache`](@ref). The clock is
advanced once for the whole step.
"""
function timestep!(integrator::ModelIntegrator, timestepper::IMEX, Δt)
    cache = integrator.state.timestepper_cache
    # type-stable split of the prognostic variables by their resolved class
    explicit_names = prognostic_names(integrator.state, cache, Val(:explicit))
    implicit_names = prognostic_names(integrator.state, cache, Val(:implicit))
    isempty(explicit_names) || timestep!(integrator, timestepper.explicit, Δt, explicit_names)
    isempty(implicit_names) || timestep!(integrator, timestepper.implicit, Δt, implicit_names)
    # advance the clock once for the entire step
    tick!(integrator.state.clock, Δt)
    return nothing
end

# Select, at compile time, the prognostic variable names whose resolved class matches `class`, by zipping
# the state's prognostic names with the IMEX cache's resolved class list (both in prognostic-variable order).
@generated function prognostic_names(::StateVariables{NF, prognames}, ::IMEXCache{classes}, ::Val{class}) where {NF, prognames, classes, class}
    selected = Symbol[name for (name, c) in zip(prognames, classes) if c === class]
    return Expr(:tuple, map(QuoteNode, selected)...)
end
