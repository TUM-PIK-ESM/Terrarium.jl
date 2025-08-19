# Base type for processes
"""
    AbstractProcess{NF}

Base type for all processes which define physics or parameterizations but are not standalone models.
"""
abstract type AbstractProcess{NF} end

variables(process::AbstractProcess) = ()

compute_auxiliary!(state, model, ::AbstractProcess) = nothing

compute_tendencies!(state, model, ::AbstractProcess) = nothing

"""
    AbstractInteraction

Base type for interactions between models or processes.
"""
abstract type AbstractInteraction end
