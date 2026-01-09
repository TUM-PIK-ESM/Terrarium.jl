abstract type AbstractHumidity end

abstract type AbstractPrecipitation end

abstract type AbstractIncomingRadiation end

abstract type AbstractAerodynamics end

"""
    $TYPEDEF

Base type for representations of the atmosphere that provide meterological state
variables such as air temperature and pressure, humidity, precipitation, incoming
solar radiation, tracer gas concentrations, wind speed, and near-surface aerodynamics.
"""
abstract type AbstractAtmosphere{
    PR<:AbstractPrecipitation,
    IR<:AbstractIncomingRadiation,
    HM<:AbstractHumidity,
    AD<:AbstractAerodynamics
} <: AbstractProcess
end
