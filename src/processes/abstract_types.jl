# Base type for processes
"""
    AbstractLandProcess{NF}

Base type for all land processes which define physics or parameterizations but are not standalone models.
"""
abstract type AbstractLandProcess{NF} end

"""
    $SIGNATURES
"""
compute_auxiliary!(idx, state, model, process::AbstractLandProcess, args...) = nothing

"""
    $SIGNATURES
"""
compute_tendencies!(idx, state, model, process::AbstractLandProcess, args...) = nothing
