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
    Humid<:AbstractHumidity,
    Precip<:AbstractPrecipitation,
    IncomingRad<:AbstractIncomingRadiation,
    Tracers<:NamedTuple{tracernames,<:Tuple{Vararg{TracerGas}}},
    Grid<:AbstractLandGrid{NF},
} <: AbstractAtmosphere{Precip, IncomingRad, Humid}
    "Spatial grid"
    grid::Grid

    "Surface-relative altitude in meters at which the atmospheric forcings are assumed to be applied"
    altitude::NF = 2*one(eltype(grid)) # Default to 2 m

    "Specific or relative humidity"
    humidity::Humid = SpecificHumidity()

    "Precipitation inputs"
    precip::Precip = TwoPhasePrecipitation()

    "Downwelling radiation inputs"
    radiation::IncomingRad = LongShortWaveRadiation()

    "Atmospheric tracer gases"
    tracers::Tracers = TracerGases(AmbientCO2())
end

variables(atmos::PrescribedAtmosphere) = (
    input(:T_air, XY(), units=u"°C", desc="Near-surface air temperature in °C"),
    input(:air_pressure, XY(), units=u"Pa", desc="Atmospheric pressure at the surface in Pa"),
    input(:windspeed, XY(), units=u"m/s", desc="Wind speed in m/s"),
    variables(atmos.humidity)...,
    variables(atmos.precip)...,
    variables(atmos.radiation)...,
    # splat all tracer variables into one tuple
    tuplejoin(map(variables, atmos.tracers)...)...,
)

"""
    air_temperature(idx, state, ::PrescribedAtmosphere)

Retrieve or compute the air temperature at the current time step.
"""
@inline air_temperature(idx, state, ::PrescribedAtmosphere) = state.air_temperature[idx]

"""
    air_pressure(idx, state, ::PrescribedAtmosphere)

Retrieve or compute the air pressure at the current time step.
"""
@inline air_pressure(idx, state, ::PrescribedAtmosphere) = state.air_pressure[idx]

"""
    windspeed(idx, state, ::PrescribedAtmosphere)

Retrieve or compute the windspeed at the current time step.
"""
@inline windspeed(idx, state, ::PrescribedAtmosphere) = state.windspeed[idx]

struct SpecificHumidity <: AbstractHumidity end

variables(::SpecificHumidity) = (
    input(:specific_humidity, XY(), units=u"kg/kg", desc="Near-surface specific humidity in kg/kg"),
)

"""
    specific_humidity(idx, state, ::PrescribedAtmosphere{PR, IR, <:SpecificHumidity})

Retrieve or compute the specific_humidity at the current time step.
"""
@inline specific_humidity(idx, state, ::PrescribedAtmosphere{PR, IR, <:SpecificHumidity}) where {PR, IR} = state.specific_humidity[idx]


struct TwoPhasePrecipitation <: AbstractPrecipitation end

variables(::TwoPhasePrecipitation) = (
    input(:rainfall, XY(), units=u"m/s", desc="Liquid precipitation (rainfall) rate"),
    input(:snowfall, XY(), units=u"m/s", desc="Frozen precipitation (snowfall) rate"),
)

"""
    rainfall(idx, state, ::PrescribedAtmosphere{<:TwoPhasePrecipitation})

Retrieve or compute the liquid precipitation (rainfall) at the current time step.
"""
@inline rainfall(idx, state, ::PrescribedAtmosphere{<:TwoPhasePrecipitation}) = state.rainfall[idx]

"""
    snowfall(idx, state, ::PrescribedAtmosphere{<:TwoPhasePrecipitation})

Retrieve or compute the frozen precipitation (snowfall) at the current time step.
"""
@inline snowfall(idx, state, ::PrescribedAtmosphere{<:TwoPhasePrecipitation}) = state.snowfall[idx]


struct LongShortWaveRadiation <: AbstractIncomingRadiation end

variables(::LongShortWaveRadiation) = (
    input(:SwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) shortwave solar radiation"),
    input(:LwIn, XY(), units=u"W/m^2", desc="Incoming (downwelling) longwave thermal radiation"),
)

"""
    shortwave_in(idx, state, ::PrescribedAtmosphere{PR, <:LongShortWaveRadiation})

Retrieve or compute the frozen precipitation (snowfall) at the current time step.
"""
shortwave_in(idx, state, ::PrescribedAtmosphere{PR, <:LongShortWaveRadiation}) where {PR} = state.SwIn[idx...]

"""
    longwave_in(idx, state, ::PrescribedAtmosphere{PR, <:LongShortWaveRadiation})

Retrieve or compute the frozen precipitation (snowfall) at the current time step.
"""
longwave_in(idx, state, ::PrescribedAtmosphere{PR, <:LongShortWaveRadiation}) where {PR} = state.LwIn[idx...]
