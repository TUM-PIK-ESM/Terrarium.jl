# Surface hydrology process types

"""
Base type for coupled surface hydrology processes.
"""
abstract type AbstractSurfaceHydrology{NF} <: AbstractCoupledProcesses{NF} end

get_evapotranspiration(hydrology::AbstractSurfaceHydrology) = hydrology.evapotranspiration

"""
Base type for canopy interception process implementations.
"""
abstract type AbstractCanopyInterception{NF} <: AbstractProcess{NF} end

"""
    canopy_water(i, j, grid, fields, ::AbstractCanopyInterception)

Compute or retrieve the current canopy water storage [kg/m^2].
"""
function canopy_water end

"""
    saturation_canopy_water(i, j, grid, fields, ::AbstractCanopyInterception)

Compute or retrieve the current canopy water saturation fraction [-].
"""
function saturation_canopy_water end

"""
    ground_precipitation(i, j, grid, fields, ::AbstractCanopyInterception)

Compute or retrieve the current rate of precipitation reaching the ground [m/s].
"""
function ground_precipitation end

"""
Base type for evapotranspiration processes.
"""
abstract type AbstractEvapotranspiration{NF} <: AbstractProcess{NF} end

"""
    surface_humidity_flux(i, j, grid, fields, ::AbstractEvapotranspiration)

Compute the surface humidity flux [m/s] at cell `i, j` based on the current state.
"""
function surface_humidity_flux end

"""
Base type for surface runoff processes.
"""
abstract type AbstractSurfaceRunoff{NF} <: AbstractProcess{NF} end

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
