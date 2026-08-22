"""
Base type for evapotranspiration processes.
"""
abstract type AbstractEvapotranspiration{NF} <: AbstractProcess{NF} end

"""
    ground_evapotranspiration_flux(i, j, grid, fields, ::AbstractEvapotranspiration)

Return the total ground evapotranspiration flux [m/s], i.e. ground evaporation + plant transpiration,
at cell `i, j` based on the current state.
"""
function ground_evapotranspiration_flux end

# Parameterizations

"""
Base type for evaporation resistance parameterizations.
"""
abstract type AbstractGroundEvaporationResistanceFactor end

"""
    ground_evaporation_resistance_factor(i, j, grid, fields, :AbstractGroundEvaporationResistanceFactor, args...)

Compute the resistance factor against ground evaporation [-] based on the current state and implementation-specific
process dependencies in `args`.
"""
function ground_evaporation_resistance_factor end
