# Surface hydrology process types

"""
Base type for canopy hydrology processes.
"""
abstract type AbstractCanopyHydrology <: AbstractProcess end

"""
Base type for evapotranspiration processes.
"""
abstract type AbstractEvapotranspiration <: AbstractProcess end

"""
Base type for surface runoff processes.
"""
abstract type AbstractRunoff <: AbstractProcess end

# Parameterizations

"""
Base type for evaporation resistance parameterizations.
"""
abstract type AbstractEvaporationResistance end
