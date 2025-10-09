# Interface for differential operators
abstract type AbstractOperator end

"""
    get_closure(op::AbstractOperator)

Returns an `AbstractClosureRelation` for the given differential operator.
Deefaults to returning `nothing` (i.e. no closure).
"""
get_closure(op::AbstractOperator)::AbstractClosureRelation = nothing

variables(op::AbstractOperator) = ()

# Interface for processes

"""
    AbstractProcess{NF}

Base type for all processes which define physics or parameterizations but are not standalone models.
"""
abstract type AbstractProcess{NF} end

variables(process::AbstractProcess) = ()

compute_auxiliary!(state, model, ::AbstractProcess) = nothing

compute_tendencies!(state, model, ::AbstractProcess) = nothing
