# Soil process types

# TODO: Think more about these process types and what methods they should have.
abstract type AbstractStratigraphy{NF} <: AbstractProcess{NF} end

abstract type AbstractSoilEnergyBalance{NF} <: AbstractProcess{NF} end

abstract type AbstractSoilHydrology{NF} <: AbstractProcess{NF} end

abstract type AbstractSoilBiogeochemistry{NF} <: AbstractProcess{NF} end
