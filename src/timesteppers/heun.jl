"""
    $TYPEDEF

Simple forward 2nd order Heun / improved Euler time stepping scheme.
"""
@kwdef struct Heun{NF} <: AbstractExplicitTimestepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0
end

Heun(::Type{NF}; kwargs...) where {NF} = Heun{NF}(; kwargs...)

default_dt(heun::Heun) = heun.Δt

is_adaptive(heun::Heun) = false

"""
    $TYPEDEF

Cache for the [`Heun`](@ref) scheme, holding copies of the prognostic state `u₀` and the predictor
tendencies `∂u∂t₀` (Heun steps in-place on `state`, so only these two are needed).
"""
struct HeunCache{NF, P, T} <: AbstractTimeStepperCache{NF}
    prognostic::P
    tendencies::T
end

# Allocate the Heun cache by deep-copying the prognostic and tendency field containers.
function initialize(::Heun, state::AbstractStateVariables)
    prognostic = map(deepcopy, state.prognostic)
    tendencies = map(deepcopy, state.tendencies)
    return HeunCache{eltype(state), typeof(prognostic), typeof(tendencies)}(prognostic, tendencies)
end

# Save the prognostic/tendency fields named in `names` into the Heun cache.
function save_cache!(cache::HeunCache, state::StateVariables, names::Tuple)
    for name in names
        copyto!(cache.prognostic[name], state.prognostic[name])
        copyto!(cache.tendencies[name], state.tendencies[name])
    end
    return nothing
end

# Restore the prognostic fields named in `names` from the Heun cache.
function restore_prognostic!(cache::HeunCache, state::StateVariables, names::Tuple)
    for name in names
        copyto!(state.prognostic[name], cache.prognostic[name])
    end
    return nothing
end

# Average the tendencies named in `names` with the saved predictor tendencies in-place:
# state.tendencies ← (state.tendencies + cache.tendencies) / 2.
function average_tendencies!(cache::HeunCache, state::StateVariables, names::Tuple)
    for name in names
        state.tendencies[name] .= (state.tendencies[name] .+ cache.tendencies[name]) ./ 2
    end
    return nothing
end

# Step only the prognostic variables named in `names` with the Heun scheme.
function timestep!(integrator::ModelIntegrator, timestepper::Heun, Δt, names::Tuple)
    (; model, state, inputs) = integrator
    grid = get_grid(model)
    # fetch Heun's cache (the whole timestepper_cache for a single Heun, or the explicit slot under IMEX)
    cache = get_cache(state.timestepper_cache, timestepper)

    # Predictor: compute tendencies ∂u∂t₀ at the current state u₀
    update_state!(state, model, inputs, compute_tendencies = true)
    # Cache u₀ and ∂u∂t₀ for the assigned variables
    save_cache!(cache, state, names)

    # Step state forward in place using ∂u∂t₀
    explicit_step!(state, grid, timestepper, Δt, names)
    # Call timestep! on model
    timestep!(state, model, timestepper, Δt)
    # Apply closure relations
    closure!(state, model)

    # Recompute tendencies ∂u∂t₁ at the predictor state (u_pred, t + Δt)
    update_state!(state, model, inputs, compute_tendencies = true)
    # Average tendencies in place: state.tendencies ← (∂u∂t₀ + ∂u∂t₁) / 2
    average_tendencies!(cache, state, names)
    # Restore the prognostic state to u₀ before the corrector explicit step
    restore_prognostic!(cache, state, names)

    # Corrector: u ← u₀ + Δt · averaged tendencies
    explicit_step!(state, grid, timestepper, Δt, names)
    # Call timestep! on model
    timestep!(state, model, timestepper, Δt)
    # Apply closure relations
    closure!(state, model)
    return nothing
end
