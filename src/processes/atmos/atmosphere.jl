abstract type AbstractSolarRadiation end

abstract type AbstractPrecipitation end

"""
Generic type representing the concentration of a particular tracer gas in the atmosphere.
"""
@kwdef struct TracerGas{name, T}
    "Tracer gas concentration in ppm"
    conc::T = nothing

    TracerGas(name::Symbol, conc) = new{name, typeof(conc)}(conc)
end

Base.nameof(::TracerGas{name}) where {name} = name

variables(::TracerGas{name}) where {name} = (
    auxiliary(name, XY(), units=u"ppm", desc="Ambient atmospheric $(name) concentration in ppm"),
)

function compute_auxiliary!(state, model, gas::TracerGas)
    conc = getproperty(state, nameof(gas))
    set!(conc, gas.conc)
end

"""
Creates a `TracerGas` for ambient CO2 with a prescribed concentration `conc`.
"""
AmbientCO2(conc=0.0) = TracerGas(:CO2, conc)
AmbientCO2(::Type{NF}) where {NF} = AmbientCO2(zero(NF))

"""
Creates a `NamedTuple` from the given tracer gas types.
"""
TracerGases(tracers::TracerGas...) = (; map(tracer -> nameof(tracer) => tracer, tracers)...)

@kwdef struct AtmosphericState{
    NF,
    tracervars,
    AirTemp,
    Humidity,
    Pressure,
    Windspeed,
    Precip<:AbstractPrecipitation,
    SolarRad<:AbstractSolarRadiation,
    Tracers<:NamedTuple{tracervars,<:Tuple{Vararg{TracerGas}}},
    Grid<:AbstractLandGrid{NF},
} <: AbstractBoundaryConditions
    "Spatial grid"
    grid::Grid

    "Surface-relative altitude in meters at which the atmospheric forcings are assumed to be applied"
    altitude::NF = 2*one(eltype(grid))
    
    "Near-surface air temperature in °C"
    T_air::AirTemp = zero(eltype(grid))

    "Near-surface specific humidity in kg/kg"
    humidity::Humidity = zero(eltype(grid))

    "Near-surface atmospheric pressure in Pa"
    pressure::Pressure = zero(eltype(grid))

    "Non-directional windspeed in m/s"
    windspeed::Windspeed = zero(eltype(grid))

    "Precipitation scheme"
    precip::Precip = TwoPhasePrecipitation(eltype(grid))

    "Downwelling solar radiation scheme"
    solar::SolarRad = TwoBandSolarRadiation(eltype(grid))

    "Atmospheric tracer gases"
    tracers::Tracers = TracerGases(AmbientCO2(zero(eltype(grid))))
end

variables(atm::AtmosphericState) = (
    auxiliary(:T_air, XY(), units=u"°C", desc="Air temperature"),
    auxiliary(:humidity, XY(), units=u"kg/kg", desc="Specific humidity"),
    auxiliary(:pressure, XY(), units=u"Pa", desc="Atmospheric pressure at the surface"),
    auxiliary(:windspeed, XY(), units=u"m/s", desc="Wind speed"),
    variables(atm.precip)...,
    variables(atm.solar)...,
    # splat all tracer variables into tuple
    tuplejoin(map(variables, atm.tracers)...)...,
)

function compute_auxiliary!(state, model, atm::AtmosphericState)
    set!(state.T_air, atm.T_air)
    set!(state.humidity, atm.humidity)
    set!(state.pressure, atm.pressure)
    set!(state.windspeed, atm.windspeed)
    compute_auxiliary!(state, model, atm.precip)
    compute_auxiliary!(state, model, atm.solar)
    for tracer in atm.tracers
        compute_auxiliary!(state, model, tracer)
    end
    return nothing
end

@kwdef struct TwoPhasePrecipitation{RF, SF} <: AbstractPrecipitation
    rainfall::RF = 0.0
    snowfall::SF = 0.0
end

TwoPhasePrecipitation(::Type{NF}) where {NF} = TwoPhasePrecipitation(zero(NF), zero(NF))

variables(::TwoPhasePrecipitation) = (
    auxiliary(:rainfall, XY(), units=u"m/s", desc="Liquid precipitation (rainfall) rate"),
    auxiliary(:snowfall, XY(), units=u"m/s", desc="Frozen precipitation (snowfall) rate"),
)

function compute_auxiliary!(state, model, precip::TwoPhasePrecipitation)
    set!(state.rainfall, precip.rainfall)
    set!(state.snowfall, precip.snowfall)
    return nothing
end

@kwdef struct TwoBandSolarRadiation{SW, LW} <: AbstractSolarRadiation
    shortwave::SW = 0.0
    longwave::LW = 0.0
end

TwoBandSolarRadiation(::Type{NF}) where {NF} = TwoBandSolarRadiation(zero(NF), zero(NF))

variables(::TwoBandSolarRadiation) = (
    auxiliary(:SwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) shortwave radiation"),
    auxiliary(:LwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) longwave radiation"),
)

function compute_auxiliary!(state, model, rad::TwoBandSolarRadiation)
    set!(state.SwIn, rad.shortwave)
    set!(state.LwIn, rad.longwave)
    return nothing
end
