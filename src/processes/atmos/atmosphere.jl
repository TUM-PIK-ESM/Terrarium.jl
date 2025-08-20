abstract type AbstractSolarRadiation end

abstract type AbstractPrecipitation end

"""
Generic type representing the concentration of a particular tracer gas in the atmosphere.
"""
@kwdef struct TracerGas{T}
    "Name of the tracer gas"
    name::Symbol

    "Tracer gas concentration in ppm"
    conc::T = nothing
end

variables(gas::TracerGas) = (
    auxiliary(gas.name, XY(), units=u"ppm", desc="Ambient atmospheric $(gas.name) concentration in ppm"),
)

function compute_auxiliary!(state, model, gas::TracerGas)
    conc = getproperty(state, gas.name)
    set!(conc, gas.conc)
end

"""
Creates a `TracerGas` for ambient CO2 with a prescribed concentration `conc`.
"""
AmbientCO2(conc=nothing) = TracerGas(:CO2, conc)

"""
Creates a `NamedTuple` from the given tracer gas types.
"""
TracerGases(tracers::TracerGas...) = (; map(tracer -> tracer.name => tracer, tracers)...)

@kwdef struct AtmosphericState{
    NF,
    tracervars,
    Grid<:AbstractLandGrid{NF},
    AirTemp,
    Humidity,
    Pressure,
    Windspeed,
    Precip<:AbstractPrecipitation,
    SolarRad<:AbstractSolarRadiation,
    Tracers<:NamedTuple{tracervars,<:Tuple{Vararg{TracerGas}}},
} <: AbstractBoundaryConditions
    "Spatial grid"
    grid::Grid

    "Surface-relative altitude in meters at which the atmospheric forcings are assumed to be applied"
    altitude::NF = 2.0
    
    "Near-surface air temperature in °C"
    T_air::AirTemp = nothing

    "Near-surface specific humidity in kg/kg"
    humidity::Humidity = nothing

    "Near-surface atmospheric pressure in Pa"
    pressure::Pressure = nothing

    "Non-directional windspeed in m/s"
    windspeed::Windspeed = nothing

    "Precipitation scheme"
    precip::Precip = TwoPhasePrecipitation()

    "Downwelling solar radiation scheme"
    solar::SolarRad = TwoBandSolarRadiation()

    "Atmospheric tracer gases"
    tracers::Tracers = TracerGases(AmbientCO2())
end

variables(atm::AtmosphericState) = (
    auxiliary(:T_air, XY(), units=u"°C", desc="Air temperature"),
    auxiliary(:humidity, XY(), units=u"kg/kg", desc="Specific humidity"),
    auxiliary(:pres, XY(), units=u"Pa", desc="Atmospheric pressure at the surface"),
    auxiliary(:wind_u, XY(), units="m/s", desc="Wind speed u-component"),
    auxiliary(:wind_v, XY(), units="m/s", desc="Wind speed v-component"),
    variables(atm.precip)...,
    variables(atm.solar)...,
    # splat all tracer variables into tuple
    tuplejoin(map(variables, atm.tracers)...)...,
)

function compute_auxiliary!(state, model, atm::AtmosphericState)
    set!(state.T_air, atm.T_air)
    set!(state.humidity, atm.humidity)
    set!(state.pres, atm.pres)
    set!(state.wind_u, atm.wind_u)
    set!(state.wind_v, atm.wind_v)
    compute_auxiliary!(state, model, atm.precip)
    compute_auxiliary!(state, model, atm.solar)
    for tracer in atm.tracers
        compute_auxiliary!(state, model, tracer)
    end
end

@kwdef struct TwoPhasePrecipitation{RF, SF} <: AbstractPrecipitation
    rainfall::RF = nothing
    snowfall::SF = nothing
end

variables(::TwoPhasePrecipitation) = (
    auxiliary(:rainfall, XY(), units=u"m/s", desc="Liquid precipitation (rainfall) rate"),
    auxiliary(:snowfall, XY(), units=u"m/s", desc="Frozen precipitation (snowfall) rate"),
)

function compute_auxiliary!(state, model, precip::TwoPhasePrecipitation)
    set!(state.rainfall, precip.rainfall)
    set!(state.snowfall, precip.snowfall)
end

@kwdef struct TwoBandSolarRadiation{SW, LW} <: AbstractSolarRadiation
    shortwave::SW = nothing
    longwave::LW = nothing
end

variables(::TwoBandSolarRadiation) = (
    auxiliary(:SwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) shortwave radiation"),
    auxiliary(:LwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) longwave radiation"),
)

function compute_auxiliary!(state, model, rad::TwoBandSolarRadiation)
    set!(state.SwIn, rad.SwIn)
    set!(state.LwIn, rad.LwIn)
end
