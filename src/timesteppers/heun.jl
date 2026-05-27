"""
    $TYPEDEF

Simple forward 2nd order Heun / improved Euler time stepping scheme.
"""
@kwdef struct Heun{NF} <: AbstractTimeStepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0
end

Heun(::Type{NF}; kwargs...) where {NF} = Heun{NF}(; kwargs...)

default_dt(heun::Heun) = heun.Δt

is_adaptive(heun::Heun) = false

# Allocate cache for the prognostic state `u₀` and the predictor tendencies `∂u∂t₀`.
# Heun steps in-place on `state`, so only these two NamedTuples of `Field`s are needed.
function initialize(::Heun, state::AbstractStateVariables)
    prognostic = map(deepcopy, state.prognostic)
    tendencies = map(deepcopy, state.tendencies)
    return (; prognostic, tendencies)
end

# Save current prognostic and tendencies into the Heun cache. Recurses into namespaces.
function save_heun_cache!(state::StateVariables)
    for name in prognostic_names(state)
        copyto!(state.cache.prognostic[name], state.prognostic[name])
        copyto!(state.cache.tendencies[name], state.tendencies[name])
    end
    fastiterate(state.namespaces) do ns
        save_heun_cache!(ns)
    end
    return nothing
end

# Restore prognostic fields from the cache. Recurses into namespaces.
function restore_heun_prognostic!(state::StateVariables)
    for name in prognostic_names(state)
        copyto!(state.prognostic[name], state.cache.prognostic[name])
    end
    fastiterate(state.namespaces) do ns
        restore_heun_prognostic!(ns)
    end
    return nothing
end

# Average current tendencies with the saved predictor tendencies in-place:
# state.tendencies ← (state.tendencies + cache.tendencies) / 2. Recurses into namespaces.
function average_heun_tendencies!(state::StateVariables)
    for name in prognostic_names(state)
        state.tendencies[name] .= (state.tendencies[name] .+ state.cache.tendencies[name]) ./ 2
    end
    fastiterate(state.namespaces) do ns
        average_heun_tendencies!(ns)
    end
    return nothing
end

function timestep!(integrator::ModelIntegrator, timestepper::Heun, Δt = default_dt(timestepper))
    (; model, state, inputs) = integrator
    grid = get_grid(model)

    # Predictor: compute tendencies ∂u∂t₀ at the current state u₀
    update_state!(state, model, inputs, compute_tendencies = true)
    # Snapshot u₀ and ∂u∂t₀
    save_heun_cache!(state)

    # Step state forward in place using ∂u∂t₀
    explicit_step!(state, grid, timestepper, Δt)
    # Call timestep! on model
    timestep!(state, model, timestepper, Δt)
    # Apply closure relations
    closure!(state, model)
 
    # Recompute tendencies ∂u∂t₁ at the predictor state (u_pred, t + Δt)
    update_state!(state, model, inputs, compute_tendencies = true)
    # Average tendencies in place: state.tendencies ← (∂u∂t₀ + ∂u∂t₁) / 2
    average_heun_tendencies!(state)
    # Restore the prognostic state to u₀ before the corrector explicit step
    restore_heun_prognostic!(state)

    # Corrector: u ← u₀ + Δt · averaged tendencies
    explicit_step!(state, grid, timestepper, Δt)
    # Call timestep! on model
    timestep!(state, model, timestepper, Δt)
    # Apply closure relations
    closure!(state, model)
    # Update clock
    tick!(state.clock, Δt)
    return nothing
end
