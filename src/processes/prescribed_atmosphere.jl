abstract type AbstractSolarRadiation end

abstract type AbstractPrecipitation end

"""
Generic type representing the concentration of a particular tracer gas in the atmosphere.
"""
@kwdef struct TracerGas{name}
    TracerGas(name::Symbol) = new{name}()
end

Base.nameof(::TracerGas{name}) where {name} = name

variables(::TracerGas{name}) where {name} = (
    input(name, XY(), units=u"ppm", desc="Ambient atmospheric $(name) concentration in ppm"),
)

"""
Creates a `TracerGas` for ambient CO2 with a prescribed concentration `conc`.
"""
AmbientCO2() = TracerGas(:CO2)

"""
Creates a `NamedTuple` from the given tracer gas types.
"""
TracerGases(tracers::TracerGas...) = (; map(tracer -> nameof(tracer) => tracer, tracers)...)

"""
"""
@kwdef struct PrescribedAtmosphere{
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
    altitude::NF = 2*one(eltype(grid)) # Default to 2 m

    "Precipitation scheme"
    precip::Precip = TwoPhasePrecipitation()

    "Downwelling solar radiation scheme"
    solar::SolarRad = TwoBandSolarRadiation()

    "Atmospheric tracer gases"
    tracers::Tracers = TracerGases(AmbientCO2())
end

variables(atm::PrescribedAtmosphere) = (
    input(:T_air, XY(), units=u"°C", desc="Near-surface air temperature in °C"),
    input(:humidity, XY(), units=u"kg/kg", desc="Near-surface specific humidity in kg/kg"),
    input(:pressure, XY(), units=u"Pa", desc="Atmospheric pressure at the surface in Pa"),
    input(:windspeed, XY(), units=u"m/s", desc="Wind speed in m/s"),
    variables(atm.precip)...,
    variables(atm.solar)...,
    # splat all tracer variables into tuple
    tuplejoin(map(variables, atm.tracers)...)...,
)

struct TwoPhasePrecipitation <: AbstractPrecipitation end

variables(::TwoPhasePrecipitation) = (
    input(:rainfall, XY(), units=u"m/s", desc="Liquid precipitation (rainfall) rate"),
    input(:snowfall, XY(), units=u"m/s", desc="Frozen precipitation (snowfall) rate"),
)

struct TwoBandSolarRadiation <: AbstractSolarRadiation end

variables(::TwoBandSolarRadiation) = (
    input(:SwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) shortwave radiation"),
    input(:LwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) longwave radiation"),
)
