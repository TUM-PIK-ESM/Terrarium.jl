# Parameterization types

"""
    $TYPEDEF

Base type for representations of near-surface atmospheric humidity. Subtypes define
which form of humidity (e.g. specific or relative humidity) is used as input.
"""
abstract type AbstractHumidity end

"""
    $TYPEDEF

Base type for representations of atmospheric precipitation. Subtypes define which
variables describe the precipitation input (e.g. rain and snow as separate fields).
"""
abstract type AbstractPrecipitation end

"""
    $TYPEDEF

Base type for representations of downwelling (incoming) radiation. Subtypes define which
spectral components are provided as inputs (e.g. split into shortwave and longwave).
"""
abstract type AbstractIncomingRadiation end

"""
    $TYPEDEF

Base type for aerodynamic parameterizations that compute the bulk drag coefficient
for turbulent heat and moisture exchange between the land surface and atmosphere.
"""
abstract type AbstractAerodynamics{NF} end

# Process types

"""
    $TYPEDEF

Base type for representations of the atmosphere that provide meterological state
variables such as air temperature and pressure, humidity, precipitation, incoming
solar radiation, tracer gas concentrations, wind speed, and near-surface aerodynamics.
"""
abstract type AbstractAtmosphere{
    NF,
    PR <: AbstractPrecipitation,
    IR <: AbstractIncomingRadiation,
    HM <: AbstractHumidity,
    AD <: AbstractAerodynamics{NF},
} <: AbstractProcess{NF}
end
