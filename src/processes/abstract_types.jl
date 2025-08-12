# Base type for processes
"""
    AbstractLandProcess{NF}

Base type for all land processes which define physics or parameterizations but are not standalone models.
"""
abstract type AbstractLandProcess{NF} end

"""
    AbstractInteraction

Base type for interactions between models or processes.
"""
abstract type AbstractInteraction end
