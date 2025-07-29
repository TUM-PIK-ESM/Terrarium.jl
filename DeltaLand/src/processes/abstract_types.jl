# Base type for processes
"""
    AbstractLandProcess

Base type for all land processes which define physics or parameterizations but are not standalone models.
"""
abstract type AbstractLandProcess end

# TODO: Think more about these process types and what methods they should have
abstract type AbstractStratigraphy <: AbstractLandProcess end
abstract type AbstractSoilEnergyBalance <: AbstractLandProcess end
abstract type AbstractSoilHydrology <: AbstractLandProcess end
abstract type AbstractSoilBiogeochemistry <: AbstractLandProcess end
