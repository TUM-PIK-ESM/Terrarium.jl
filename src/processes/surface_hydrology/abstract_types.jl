# Surface hydrology process types

"""
Base type for canopy hydrology processes.
"""
abstract type AbstractCanopyHydrology <: AbstractProcess end

"""
    canopy_water(i, j, state, grid, ::AbstractCanopyHydrology)

Compute or retrieve the current canopy water storage [kg/m^2].
"""
function canopy_water end

"""
    saturation_canopy_water(i, j, state, grid, ::AbstractCanopyHydrology)

Compute or retrieve the current canopy water saturation fraction [-].
"""
function saturation_canopy_water end

"""
    ground_precipitation(i, j, state, grid, ::AbstractCanopyHydrology)

Compute or retrieve the current rate of precipitation reaching the ground [m/s].
"""
function ground_precipitation end

"""
Base type for evapotranspiration processes.
"""
abstract type AbstractEvapotranspiration <: AbstractProcess end

"""
    surface_humidity_flux(i, j, state, grid, ::AbstractEvapotranspiration)

Compute the surface humidity flux [m/s] at cell `i, j` based on the current state.
"""
function surface_humidity_flux end

"""
Base type for surface runoff processes.
"""
abstract type AbstractSurfaceRunoff <: AbstractProcess end

"""
Base type for coupled surface hydrology processes.
"""
abstract type AbstractSurfaceHydrology <: AbstractProcess end

# Parameterizations

"""
Base type for evaporation resistance parameterizations.
"""
abstract type AbstractGroundEvaporationResistanceFactor end

"""
    ground_evaporation_resistance_factor(i, j, state, grid, :AbstractGroundEvaporationResistanceFactor, args...)

Compute the resistance factor against ground evaporation [-] based on the current state and implementation-specific
process dependencies in `args`.
"""
function ground_evaporation_resistance_factor end
