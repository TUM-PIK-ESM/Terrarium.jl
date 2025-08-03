# Soil process types

# TODO: Think more about these process types and what methods they should have.
# Also maybe move into respective subfolders.
abstract type AbstractStratigraphy <: AbstractLandProcess end

abstract type AbstractSoilEnergyBalance <: AbstractLandProcess end

abstract type AbstractSoilHydrology <: AbstractLandProcess end

abstract type AbstractSoilBiogeochemistry <: AbstractLandProcess end
