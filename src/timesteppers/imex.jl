"""
    $TYPEDEF

Implicit-explicit (IMEX) time stepper that integrates each prognostic variable with one of two
sub-steppers depending on its timestepper class (see [`prognostic`](@ref) and the `timestepper_classes`
keyword of [`StateVariables`](@ref)): variables of class `:explicit` are stepped by `explicit`, and those
of class `:implicit` are stepped by `implicit`.

Properties:
$(TYPEDFIELDS)
"""
struct IMEX{
        NF,
        E <: AbstractExplicitTimestepper{NF},
        I <: AbstractImplicitTimestepper{NF},
    } <: AbstractTimeStepper{NF}
    "Sub-stepper for prognostic variables of class `:explicit`"
    explicit::E

    "Sub-stepper for prognostic variables of class `:implicit`"
    implicit::I
end

"""
    IMEX(; explicit::AbstractExplicitTimestepper, implicit::AbstractImplicitTimestepper)
    IMEX(explicit, implicit)

Construct an [`IMEX`](@ref) time stepper from an `explicit` and an `implicit` sub-stepper. Both must
share the same numerical (float) type.
"""
IMEX(; explicit, implicit) = IMEX(explicit, implicit)

default_dt(imex::IMEX) = default_dt(imex.explicit)

is_adaptive(imex::IMEX) = is_adaptive(imex.explicit) || is_adaptive(imex.implicit)

# An IMEX cache holds the two timestepper caches, one per sub-stepper.
initialize_timestepper_cache(imex::IMEX, state) = (; explicit = initialize(imex.explicit, state), implicit = initialize(imex.implicit, state))

# Select a sub-stepper's cache from an IMEX cache (the bare-cache fallbacks live in abstract_timestepper.jl).
explicit_cache(cache::NamedTuple{(:explicit, :implicit)}) = cache.explicit
implicit_cache(cache::NamedTuple{(:explicit, :implicit)}) = cache.implicit

"""
    timestep!(integrator::ModelIntegrator, timestepper::IMEX, Δt)

Advance the model forward by one timestep of size `Δt` using an [`IMEX`](@ref) timestepper, integrating
each prognostic variable with the sub-stepper matching its timestepper class (`:explicit`/`:implicit`).
The prognostic variables are grouped by class and each group is forwarded to the corresponding
`timestep!(integrator, sub_timestepper, Δt, names)` method. The clock is advanced once for the whole step.
"""
function timestep!(integrator::ModelIntegrator, timestepper::IMEX, Δt)
    # step the explicit-class prognostic variables with the explicit sub-stepper
    explicit_names = prognostic_names(integrator.state, :explicit)
    isempty(explicit_names) || timestep!(integrator, timestepper.explicit, Δt, explicit_names)
    # step the implicit-class prognostic variables with the implicit sub-stepper
    implicit_names = prognostic_names(integrator.state, :implicit)
    isempty(implicit_names) || timestep!(integrator, timestepper.implicit, Δt, implicit_names)
    # advance the clock once for the entire step
    tick!(integrator.state.clock, Δt)
    return nothing
end
