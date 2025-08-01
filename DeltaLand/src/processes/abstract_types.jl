# Base type for processes
"""
    AbstractLandProcess

Base type for all land processes which define physics or parameterizations but are not standalone models.
"""
abstract type AbstractLandProcess end

"""
    $SIGNATURES
"""
update_state!(idx, state, model::AbstractModel, process::AbstractLandProcess, args...) = nothing

"""
    $SIGNATURES
"""
compute_tendencies!(idx, state, model::AbstractModel, process::AbstractLandProcess, args...) = nothing

# Soil processes

# TODO: Think more about these process types and what methods they should have.
# Also maybe move into respective subfolders.
abstract type AbstractStratigraphy <: AbstractLandProcess end

abstract type AbstractSoilEnergyBalance <: AbstractLandProcess end

abstract type AbstractSoilHydrology <: AbstractLandProcess end

abstract type AbstractSoilBiogeochemistry <: AbstractLandProcess end

# Vegetation processes

abstract type AbstractPhotosynthesis end

abstract type AbstractStomatalConductance end

abstract type AbstractVegetationCarbonDynamics end
