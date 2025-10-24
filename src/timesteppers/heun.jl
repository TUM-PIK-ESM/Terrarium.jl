"""
    $TYPEDEF

Simple forward 2nd order Heun / improved Euler time stepping scheme.
"""
@kwdef struct Heun{NF,C} <: AbstractTimeStepper{NF}
    "Initial timestep size in seconds"
    Δt::NF = 300.0

    "Cache for intermediate results"
    cache::C
end

"""
    Heun(state::StateVariables; kwargs...)

Create a `Heun` timestepper for the given state variables.
"""
Heun(state::StateVariables; kwargs...) = Heun{eltype(state), typeof(state)}(; cache=deepcopy(state), kwargs...)

default_dt(heun::Heun) = heun.Δt

is_adaptive(heun::Heun) = false

function average_tendencies!(state::StateVariables{prognames, tendnames}, cache::StateVariables{prognames, tendnames}) where {prognames, tendnames}
    for tendname in tendnames
        state.tendencies[tendname] .= (state.tendencies[tendname] + cache.tendencies[tendname]) / 2
    end
end

function timestep!(state, model::AbstractModel, timestepper::Heun, Δt = default_dt(timestepper))
    
    cache = timestepper.cache
    
    # trial step 
    compute_auxiliary!(state, model) 
    compute_tendency!(state, model) 

    copy!(cache, state)

    # improved step 
    explicit_step(cache, model, timestepper, Δt) # for 2nd order Euler you don't do a step in the middle   
    compute_auxiliary!(cache, model) 
    compute_tendency!(cache, model) 

    # improved Euler step call that steps `state` forward but averages `state.tendencies` 
    average_tendencies!(state, cache)
    explicit_step(state, model, timestepper, Δt) 
    
    # Apply inverse closure relations
    for closure in state.closures
        invclosure!(state, model, closure)
    end
    # Update clock
    tick!(state.clock, Δt)
    return nothing
end
