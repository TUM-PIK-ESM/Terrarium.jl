# Soil process types
# TODO: Think more about these process types and what methods they should have.

abstract type AbstractStratigraphy{NF} <: AbstractProcess end

abstract type AbstractSoilEnergyBalance{NF} <: AbstractProcess end

abstract type AbstractSoilHydrology{NF} <: AbstractProcess end

abstract type AbstractSoilBiogeochemistry{NF} <: AbstractProcess end
