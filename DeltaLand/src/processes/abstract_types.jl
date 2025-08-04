# Base type for processes
"""
    AbstractLandProcess

Base type for all land processes which define physics or parameterizations but are not standalone models.
"""
abstract type AbstractLandProcess end

"""
    $SIGNATURES
"""
compute_auxiliary!(idx, state, model::AbstractModel, process::AbstractLandProcess, args...) = nothing

"""
    $SIGNATURES
"""
compute_tendencies!(idx, state, model::AbstractModel, process::AbstractLandProcess, args...) = nothing
