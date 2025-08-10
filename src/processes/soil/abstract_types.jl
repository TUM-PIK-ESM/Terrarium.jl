# Soil process types

# TODO: Think more about these process types and what methods they should have.
# Also maybe move into respective subfolders.
abstract type AbstractStratigraphy{NF} <: AbstractLandProcess{NF} end

abstract type AbstractSoilEnergyBalance{NF} <: AbstractLandProcess{NF} end

abstract type AbstractSoilHydrology{NF} <: AbstractLandProcess{NF} end

abstract type AbstractSoilBiogeochemistry{NF} <: AbstractLandProcess{NF} end
