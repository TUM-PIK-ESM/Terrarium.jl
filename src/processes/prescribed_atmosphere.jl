abstract type AbstractPrecipitation end

abstract type AbstractIncomingRadiation end

abstract type AbstractAtmosphere{
    PR<:AbstractPrecipitation,
    IR<:AbstractIncomingRadiation
} <: AbstractBoundaryConditions
end

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
Creates a `TracerGas` for ambient CO2 with concentration prescribed by an input variable with
the given name.
"""
AmbientCO2(name=:CO2) = TracerGas(name)

"""
Creates a `NamedTuple` from the given tracer gas types.
"""
TracerGases(tracers::TracerGas...) = (; map(tracer -> nameof(tracer) => tracer, tracers)...)

"""
    $TYPEDEF

Represents prescribed atmospheric conditions given by the following input variables:
    - Air temperature
    - Humidity
    - Atmospheric pressure
    - Windspeed
    - Precipitation
    - Solar radiation
    - Zero or more tracer gases (defaults to CO2 only)

Precpitation and solar radiation are specified according to specialized subtypes which dictate
the form of the input data; for precipitation, this defaults to `TwoPhasePrecipitation`, i.e.
rain- and snowfall given as separate inputs, while for solar radiation, the default is
`LongShortWaveRadiation` which partitions downwelling radiation into the common short- and long
wave lengths representing solar and thermal (infrared) radiation.
"""
@kwdef struct PrescribedAtmosphere{
    NF,
    tracernames,
    AirTemp,
    Humidity,
    Pressure,
    Windspeed,
    Precip<:AbstractPrecipitation,
    IncomingRad<:AbstractIncomingRadiation,
    Tracers<:NamedTuple{tracernames,<:Tuple{Vararg{TracerGas}}},
    Grid<:AbstractLandGrid{NF},
} <: AbstractAtmosphere{Precip, IncomingRad}
    "Spatial grid"
    grid::Grid

    "Surface-relative altitude in meters at which the atmospheric forcings are assumed to be applied"
    altitude::NF = 2*one(eltype(grid)) # Default to 2 m

    "Precipitation inputs"
    precip::Precip = TwoPhasePrecipitation()

    "Downwelling radiation inputs"
    radiation::IncomingRad = LongShortWaveRadiation()

    "Atmospheric tracer gases"
    tracers::Tracers = TracerGases(AmbientCO2())
end

variables(atm::PrescribedAtmosphere) = (
    input(:T_air, XY(), units=u"°C", desc="Near-surface air temperature in °C"),
    input(:humidity, XY(), units=u"kg/kg", desc="Near-surface specific humidity in kg/kg"),
    input(:pressure, XY(), units=u"Pa", desc="Atmospheric pressure at the surface in Pa"),
    input(:windspeed, XY(), units=u"m/s", desc="Wind speed in m/s"),
    variables(atm.precip)...,
    variables(atm.radiation)...,
    # splat all tracer variables into one tuple
    tuplejoin(map(variables, atm.tracers)...)...,
)

struct TwoPhasePrecipitation <: AbstractPrecipitation end

variables(::TwoPhasePrecipitation) = (
    input(:rainfall, XY(), units=u"m/s", desc="Liquid precipitation (rainfall) rate"),
    input(:snowfall, XY(), units=u"m/s", desc="Frozen precipitation (snowfall) rate"),
)

struct LongShortWaveRadiation <: AbstractSolarRadiation end

variables(::LongShortWaveRadiation) = (
    input(:SwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) shortwave solar radiation"),
    input(:LwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) longwave thermal radiation"),
)
