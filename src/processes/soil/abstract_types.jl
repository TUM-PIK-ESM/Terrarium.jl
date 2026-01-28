# Soil process types

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

Base type for mineral soil texture parameterizations.
"""
abstract type AbstractSoilTexture{NF} end

"""
    $TYPEDEF

Base type for parameterizations of soil porosity.
"""
abstract type AbstractSoilPorosity{NF} end

"""
    mineral_porosity(::AbstractSoilPorosity, texture::SoilTexture)

Compute or retrieve the natural porosity of the mineral soil constitutents, i.e.
excluding organic material.
"""
function mineral_porosity end

"""
    organic_porosity(::AbstractSoilPorosity, texture::SoilTexture)

Compute or retrieve the natural porosity of the organic soil constitutents, i.e.
excluding mineral material.
"""
function organic_porosity end

"""
    $TYPEDEF

Base type for soil stratigraphy parameterizations.
"""
abstract type AbstractStratigraphy{NF} end

"""
    soil_texture(i, j, k, grid, state, ::AbstractStratigraphy, args...)

Return the texture of the soil at index `i, j, k` for the given stratigraphy parameterization.
"""
function soil_texture end

"""
    soil_matrix(i, j, k, grid, state, ::AbstractStratigraphy, args...)

Return the solid matrix of the soil at index `i, j, k` for the given stratigraphy parameterization.
"""
function soil_matrix end

"""
    soil_volume(i, j, k, grid, state, ::AbstractStratigraphy, args...)

Return a description of the full material composition of the soil volume at index `i, j, k` for the
given stratigraphy parameterization.
"""
function soil_volume end

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
