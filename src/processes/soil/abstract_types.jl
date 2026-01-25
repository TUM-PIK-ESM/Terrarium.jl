# Soil process types

abstract type AbstractStratigraphy{NF} <: AbstractProcess end

abstract type AbstractSoilEnergyBalance{NF} <: AbstractProcess end

"""
    energy_tendency(i, j, k, grid, ::SoilEnergyBalance, args...)

Compute the internal energy tendency `∂U∂t` at index `i, j, k`.
"""
function energy_tendency end

"""
    thermal_conductivity(i, j, k, grid, ::SoilEnergyBalance, args...)

Compute the thermal conductivity at index `i, j, k`.
"""
function thermal_conductivity end

abstract type AbstractSoilHydrology{NF} <: AbstractProcess end

abstract type AbstractSoilBiogeochemistry{NF} <: AbstractProcess end

# Parameterization types

"""
    $TYPEDEF

Base type for formulations of the heat transfer operator.
"""
abstract type AbstractHeatOperator end

"""
    $TYPEDEF

Base type for closure relations between energy and temperature in soil volumes.
"""
abstract type AbstractSoilEnergyClosure <: AbstractClosureRelation end

"""
    $TYPEDEF

Base type for closure relations between water saturation and potential in soil volumes.
"""
abstract type AbstractSoilWaterClosure <: AbstractClosureRelation end
