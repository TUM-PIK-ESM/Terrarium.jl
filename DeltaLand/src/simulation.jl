abstract type AbstractStateVariables end

abstract type AbstractSimulation end

struct Simulation{Model<:AbstractModel, StateVars<:AbstractStateVariables} <: AbstractSimulation
    model::Model
    
    state::StateVars

    clock::Clock
end

function timestep!(sim::Simulation, dt=nothing)
    timestepper = get_time_stepping(sim.model)
    dt = isnothing(dt) ? get_dt(timestepper) : dt
    timestep!(sim.state, sim.model, timestepper, dt)
end

# temporary solution for holding state variables
@kwdef struct StateVariables{prognames,tendnames,auxnames,ProgVars,TendVars,AuxVars} <: AbstractStateVariables
    prognostic::NamedTuple{prognames,ProgVars}
    tendencies::NamedTuple{tendnames,TendVars}
    auxiliary::NamedTuple{auxnames,AuxVars}
end

Base.propertynames(vars::StateVariables{prognames,tendnames,auxnames}) where {prognames,tendnames,auxnames} = (
    prognames...,
    tendnames...,
    auxnames...,
    fieldnames(vars)...,
)

function Base.getproperty(vars::StateVariables, name::Symbol)
    # forward getproperty calls to variable groups
    if name ∈ getfield(vars, :prognostic)
        return getproperty(getfield(vars, :prognostic), name)
    elseif name ∈ getfield(vars, :tendencies)
        return getproperty(getfield(vars, :tendencies), name)
    elseif name ∈ getfield(vars, :auxiliary)
        return getproperty(getfield(vars, :auxiliary), name)
    else
        return getfield(vars, name)
    end
end