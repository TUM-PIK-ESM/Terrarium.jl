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

# also allow dispatch on nothing
compute_auxiliary!(state, model, ::Nothing) = nothing
compute_tendencies!(state, model, ::Nothing) = nothing

# Interface for differential operators
abstract type AbstractOperator end

variables(op::AbstractOperator) = ()

"""
    get_closure(op::AbstractOperator)

Returns an `AbstractClosureRelation` for the given differential operator.
Deefaults to returning `nothing` (i.e. no closure).
"""
get_closure(op::AbstractOperator) = nothing

closure!(state, model::AbstractModel, ::Nothing) = nothing

invclosure!(state, model::AbstractModel, ::Nothing) = nothing
