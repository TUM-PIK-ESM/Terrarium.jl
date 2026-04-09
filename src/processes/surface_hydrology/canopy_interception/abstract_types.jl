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
    rainfall_ground(i, j, grid, fields, ::AbstractCanopyInterception)

Compute or retrieve the current rate of precipitation reaching the ground [m/s].
"""
function rainfall_ground end
