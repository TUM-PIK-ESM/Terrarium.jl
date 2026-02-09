# Parameterization types

abstract type AbstractHumidity end

abstract type AbstractPrecipitation end

abstract type AbstractIncomingRadiation end

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
