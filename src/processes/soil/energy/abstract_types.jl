"""
    $TYPEDEF

Base type for soil energy balance process implementations. Subtypes should define
state variables for soil `temperature`, `internal_energy`, and any other relevant
thermal properties or state variables. Soil energy balances evolve the soil internal
energy and therefore subtype [`AbstractEnergyBalance`](@ref).
"""
abstract type AbstractSoilEnergyBalance{NF} <: AbstractEnergyBalance{NF} end

# Parameterizations

"""
Base type for bulk weighting/mixing schemes that calculate weighted mixture of material properties
such as conductivities or densities.
"""
abstract type AbstractBulkWeighting end
