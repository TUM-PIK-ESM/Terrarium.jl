# Interface for processes

"""
    AbstractProcess

Base type for all processes which define physics or parameterizations but are not standalone models.
"""
abstract type AbstractProcess end

variables(process::AbstractProcess) = ()

initialize!(state, model, process::AbstractProcess) = compute_auxiliary!(state, model, process)

compute_auxiliary!(state, model, ::AbstractProcess) = nothing

compute_tendencies!(state, model, ::AbstractProcess) = nothing

# Also allow dispatch on nothing
compute_auxiliary!(state, model, ::Nothing) = nothing
compute_tendencies!(state, model, ::Nothing) = nothing

# Allow dispatch on nothing types for processes
closure!(state, model::AbstractModel, ::Nothing) = nothing
invclosure!(state, model::AbstractModel, ::Nothing) = nothing
