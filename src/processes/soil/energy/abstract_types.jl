"""
    $TYPEDEF

Base type for soil energy balance process implementations. Subtypes should define
state variables for soil `temperature`, `internal_energy`, and any other relevant
thermal properties or state variables.
"""
abstract type AbstractSoilEnergyBalance{NF} <: AbstractProcess{NF} end

# Method interface

"""
    get_thermal_properties(energy::AbstractSoilEnergyBalance)

Return the thermal properites associated with the given soil energy balance process.
"""
function get_thermal_properties end

"""
    compute_energy_tendency(i, j, k, grid, ::AbstractSoilEnergyBalance, args...)

Compute the internal energy tendency `∂U∂t` at index `i, j, k`.
"""
function compute_energy_tendency end

"""
    compute_thermal_conductivity(i, j, k, grid, ::AbstractSoilEnergyBalance, args...)

Compute the thermal conductivity at index `i, j, k`.
"""
function compute_thermal_conductivity end

# Process methods

"""
    initialize!(state, grid, energy::AbstractSoilEnergyBalance, soil::AbstractSoil, constants::PhysicalConstants)

Initialize energy state from the current `state` on `grid` for the given `soil` configuration. Unless otherwise stated,
it should generally be assumed that the `temperature` is already initialized prior to this method being called.
"""
initialize!(state, grid, energy::AbstractSoilEnergyBalance, soil::AbstractSoil, constants::PhysicalConstants) = nothing

"""
    compute_auxiliary!(state, grid, energy::AbstractSoilEnergyBalance, soil::AbstractSoil, constants::PhysicalConstants)

Compute energy auxiliaries from the current `state` on `grid` for the given `soil` configuration.
"""
compute_auxiliary!(state, grid, energy::AbstractSoilEnergyBalance, soil::AbstractSoil, constants::PhysicalConstants) = nothing

"""
    compute_tendencies!(state, grid, energy::AbstractSoilEnergyBalance, soil::AbstractSoil, constants::PhysicalConstants)

Compute energy tendencies from the current `state` on `grid` for the given `soil` configuration.
"""
compute_tendencies!(state, grid, energy::AbstractSoilEnergyBalance, soil::AbstractSoil, constants::PhysicalConstants) = nothing

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
