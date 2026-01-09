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
struct PrescribedAtmosphere{
    NF,
    tracernames,
    Humid<:AbstractHumidity,
    Precip<:AbstractPrecipitation,
    IncomingRad<:AbstractIncomingRadiation,
    Tracers<:NamedTuple{tracernames,<:Tuple{Vararg{TracerGas}}}
} <: AbstractAtmosphere{Precip, IncomingRad, Humid}
    "Surface-relative altitude in meters at which the atmospheric forcings are assumed to be applied"
    altitude::NF

    "Specific or relative humidity"
    humidity::Humid

    "Precipitation inputs"
    precip::Precip

    "Downwelling radiation inputs"
    radiation::IncomingRad

    "Aerodynamic resistances and drag coefficients"
    aerodynamics::Aerodynamics

    "Atmospheric tracer gases"
    tracers::Tracers
end

function PrescribedAtmosphere(
    ::Type{NF};
    altitude::NF = NF(10), # Default to 10 m
    humidity::AbstractHumidity = SpecificHumidity(),
    precip::AbstractPrecipitation = TwoPhasePrecipitation(),
    radiation::AbstractIncomingRadiation = LongShortWaveRadiation(),
    tracers::NamedTuple = TracerGases(AmbientCO2()),
) where {NF}
    return PrescribedAtmosphere(altitude, humidity, precip, radiation, tracers)
end

variables(atmos::PrescribedAtmosphere) = (
    input(:air_temperature, XY(), units=u"°C", desc="Near-surface air temperature in °C"),
    input(:air_pressure, XY(), units=u"Pa", desc="Atmospheric pressure at the surface in Pa"),
    input(:windspeed, XY(), units=u"m/s", desc="Wind speed in m/s"),
    variables(atmos.humidity)...,
    variables(atmos.precip)...,
    variables(atmos.radiation)...,
    # splat all tracer variables into one tuple
    tuplejoin(map(variables, atmos.tracers)...)...,
)

@inline compute_auxiliary!(state, model, atmos::AbstractAtmosphere) = nothing

@inline compute_tendencies!(state, model, atmos::AbstractAtmosphere) = nothing

"""
    air_temperature(i, j, state, ::AbstractAtmosphere)

Retrieve or compute the air temperature at the current time step.
"""
@inline air_temperature(i, j, state, ::AbstractAtmosphere) = state.air_temperature[i, j]

"""
    air_pressure(i, j, state, ::AbstractAtmosphere)

Retrieve or compute the air pressure at the current time step.
"""
@inline air_pressure(i, j, state, ::AbstractAtmosphere) = state.air_pressure[i, j]

"""
    windspeed(i, j, state, ::PrescribedAtmosphere)

Retrieve or compute the windspeed at the current time step.
"""
@inline windspeed(i, j, state, ::AbstractAtmosphere) = state.windspeed[i, j]

struct SpecificHumidity <: AbstractHumidity end

variables(::SpecificHumidity) = (
    input(:specific_humidity, XY(), units=u"kg/kg", desc="Near-surface specific humidity in kg/kg"),
)

"""
    specific_humidity(i, j, state, ::PrescribedAtmosphere{PR, IR, <:SpecificHumidity})

Retrieve or compute the specific_humidity at the current time step.
"""
@inline specific_humidity(i, j, state, ::AbstractAtmosphere{PR, IR, <:SpecificHumidity}) where {PR, IR} = state.specific_humidity[i, j]

"""
    $SIGNATURES

Computes the specific humidity (vapor pressure) deficit over a surface at temperature `Ts` from the current atmospheric state.
"""
@inline function compute_humidity_vpd(i, j, atmos::AbstractAtmosphere, c::PhysicalConstants, Ts = nothing)
    let Δe = compute_vpd(i, j, atmos, c, Ts),
        p = air_pressure(i, j, state, atmos);
        Δq = vapor_pressure_to_specific_humidity(Δe, p, c.ε)
        return Δq
    end
end

"""
    $SIGNATURES

Computes the vapor pressure deficit over a surface at temperature `Ts` from the current atmospheric state.
"""
@inline function compute_vpd(i, j, atmos::AbstractAtmosphere, c::PhysicalConstants, Ts = nothing) where NF
    @inbounds let Tair = air_temperature(i, j, state, atmos),
                  q_air = specific_humidity(i, j, state, atmos),
                  pres = air_pressure(i, j, state, atmos);
        return compute_vpd(c, pres, q_air, Tair, Ts)
    end
end

struct TwoPhasePrecipitation <: AbstractPrecipitation end

variables(::TwoPhasePrecipitation) = (
    input(:rainfall, XY(), units=u"m/s", desc="Liquid precipitation (rainfall) rate"),
    input(:snowfall, XY(), units=u"m/s", desc="Frozen precipitation (snowfall) rate"),
)

"""
    rainfall(i, j, state, ::AbstractAtmosphere{<:TwoPhasePrecipitation})

Retrieve or compute the liquid precipitation (rainfall) at the current time step.
"""
@inline rainfall(i, j, state, ::AbstractAtmosphere{<:TwoPhasePrecipitation}) = state.rainfall[i, j]

"""
    snowfall(i, j, state, ::AbstractAtmosphere{<:TwoPhasePrecipitation})

Retrieve or compute the frozen precipitation (snowfall) at the current time step.
"""
@inline snowfall(i, j, state, ::AbstractAtmosphere{<:TwoPhasePrecipitation}) = state.snowfall[i, j]

struct LongShortWaveRadiation <: AbstractIncomingRadiation end

variables(::LongShortWaveRadiation) = (
    input(:surface_shortwave_down, XY(), units=u"W/m^2", desc="Incoming (downwelling) shortwave solar radiation"),
    input(:surface_longwave_down, XY(), units=u"W/m^2", desc="Incoming (downwelling) longwave thermal radiation"),
)

"""
    shortwave_in(i, j, state, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation})

Retrieve or compute the incoming/downwelling shortwave radiation at the current time step.
"""
shortwave_in(i, j, state, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation}) where {PR} = state.surface_shortwave_down[i, j]

"""
    longwave_in(i, j, state, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation})

Retrieve or compute the incoming/downwelling longwave radiation at the current time step.
"""
longwave_in(i, j, state, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation}) where {PR} = state.surface_longwave_down[i, j]
