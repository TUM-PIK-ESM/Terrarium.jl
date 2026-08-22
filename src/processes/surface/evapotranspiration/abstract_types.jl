"""
Base type for evapotranspiration processes.
"""
abstract type AbstractEvapotranspiration{NF} <: AbstractProcess{NF} end

"""
    surface_evaporation_flux(i, j, grid, fields, ::AbstractEvapotranspiration)

Compute the surface evaporation flux [m/s] at cell `i, j` based on the current state.
"""
function surface_evaporation_flux end

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
\
