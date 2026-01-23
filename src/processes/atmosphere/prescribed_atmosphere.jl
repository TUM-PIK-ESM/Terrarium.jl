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
    Humidity<:AbstractHumidity,
    Precip<:AbstractPrecipitation,
    IncomingRad<:AbstractIncomingRadiation,
    Aerodynamics<:AbstractAerodynamics,
    Tracers<:NamedTuple{tracernames,<:Tuple{Vararg{TracerGas}}}
} <: AbstractAtmosphere{Precip, IncomingRad, Humidity, Aerodynamics}
    "Surface-relative altitude in meters at which the atmospheric forcings are assumed to be applied"
    altitude::NF

    "Minimum windspeed"
    min_windspeed::NF

    "Specific or relative humidity"
    humidity::Humidity

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
    min_windspeed::NF = NF(0.01), # Default to 0.01 m/s
    humidity::AbstractHumidity = SpecificHumidity(),
    precip::AbstractPrecipitation = TwoPhasePrecipitation(),
    radiation::AbstractIncomingRadiation = LongShortWaveRadiation(),
    aerodynamics::AbstractAerodynamics = ConstantAerodynamics(),
    tracers::NamedTuple = TracerGases(AmbientCO2()),
) where {NF}
    return PrescribedAtmosphere(altitude, min_windspeed, humidity, precip, radiation, aerodynamics, tracers)
end

variables(atmos::PrescribedAtmosphere{NF}) where {NF} = (
    input(:air_temperature, XY(), default=NF(10), units=u"°C", desc="Near-surface air temperature in °C"),
    input(:air_pressure, XY(), default=NF(101_325), units=u"Pa", desc="Atmospheric pressure at the surface in Pa"),
    input(:windspeed, XY(), default=NF(0.1), units=u"m/s", desc="Wind speed in m/s"),
    variables(atmos.humidity)...,
    variables(atmos.precip)...,
    variables(atmos.radiation)...,
    variables(atmos.aerodynamics)...,
    # splat all tracer variables into one tuple
    tuplejoin(map(variables, atmos.tracers)...)...,
)

@inline compute_auxiliary!(state, model, atmos::PrescribedAtmosphere) = nothing

@inline compute_tendencies!(state, model, atmos::PrescribedAtmosphere) = nothing

"""
    aerodynamic_resistance(i, j, state, grid, atmos::PrescribedAtmosphere)

Compute the aerodynamic resistance (inverse conductance) at grid cell `i, j`.
"""
@inline function aerodynamic_resistance(i, j, state, grid, atmos::PrescribedAtmosphere)
    let C = drag_coefficient(i, j, state, grid, atmos.aerodynamics),
        Vₐ = max(windspeed(i, j, state, grid, atmos), 1e-6); # clip windspeed to small value
        rₐ = 1 / (C * Vₐ)
        return rₐ
    end
end

"""
    air_temperature(i, j, state, grid, ::PrescribedAtmosphere)

Retrieve or compute the air temperature at the current time step.
"""
@propagate_inbounds air_temperature(i, j, state, grid, ::PrescribedAtmosphere) = state.air_temperature[i, j]

"""
    air_pressure(i, j, state, grid, ::PrescribedAtmosphere)

Retrieve or compute the air pressure at the current time step.
"""
@propagate_inbounds air_pressure(i, j, state, grid, ::PrescribedAtmosphere) = state.air_pressure[i, j]

"""
    windspeed(i, j, state, ::PrescribedAtmosphere)

Retrieve or compute the windspeed at the current time step.
"""
@propagate_inbounds windspeed(i, j, state, grid, atmos::PrescribedAtmosphere) = max(state.windspeed[i, j], atmos.min_windspeed)

struct SpecificHumidity <: AbstractHumidity end

variables(::SpecificHumidity) = (
    input(:specific_humidity, XY(), units=u"kg/kg", desc="Near-surface specific humidity in kg/kg"),
)

"""
    specific_humidity(i, j, state, ::PrescribedAtmosphere{PR, IR, <:SpecificHumidity})

Retrieve or compute the specific_humidity at the current time step.
"""
@propagate_inbounds specific_humidity(i, j, state, grid, ::PrescribedAtmosphere{PR, IR, <:SpecificHumidity}) where {PR, IR} = state.specific_humidity[i, j]

"""
    $TYPEDSIGNATURES

Computes the specific humidity (vapor pressure) deficit over a surface at temperature `Ts` from the current atmospheric state.
"""
@propagate_inbounds function compute_humidity_vpd(i, j, state, grid, atmos::AbstractAtmosphere, c::PhysicalConstants, Ts = nothing)
    let Δe = compute_vpd(i, j, state, grid, atmos, c, Ts),
        p = air_pressure(i, j, state, grid, atmos);
        Δq = vapor_pressure_to_specific_humidity(Δe, p, c.ε)
        return Δq
    end
end

"""
    $TYPEDSIGNATURES

Computes the vapor pressure deficit over a surface at temperature `Ts` from the current atmospheric state.
"""
@propagate_inbounds function compute_vpd(i, j, state, grid, atmos::AbstractAtmosphere, c::PhysicalConstants, Ts = nothing)
    Tair = air_temperature(i, j, state, grid, atmos)
    q_air = specific_humidity(i, j, state, grid, atmos)
    pres = air_pressure(i, j, state, grid, atmos)
    Ts = isnothing(Ts) ? Tair : Ts
    return compute_vpd(c, pres, q_air, Ts)
end

struct TwoPhasePrecipitation <: AbstractPrecipitation end

variables(::TwoPhasePrecipitation) = (
    input(:rainfall, XY(), units=u"m/s", desc="Liquid precipitation (rainfall) rate"),
    input(:snowfall, XY(), units=u"m/s", desc="Frozen precipitation (snowfall) rate"),
)

"""
    rainfall(i, j, state, grid, ::AbstractAtmosphere{<:TwoPhasePrecipitation})

Retrieve or compute the liquid precipitation (rainfall) at the current time step.
"""
@inline rainfall(i, j, state, grid, ::AbstractAtmosphere{<:TwoPhasePrecipitation}) = state.rainfall[i, j]

"""
    snowfall(i, j, state, grid, ::AbstractAtmosphere{<:TwoPhasePrecipitation})

Retrieve or compute the frozen precipitation (snowfall) at the current time step.
"""
@inline snowfall(i, j, state, grid, ::AbstractAtmosphere{<:TwoPhasePrecipitation}) = state.snowfall[i, j]

struct LongShortWaveRadiation <: AbstractIncomingRadiation end

variables(::LongShortWaveRadiation) = (
    input(:surface_shortwave_down, XY(), units=u"W/m^2", desc="Incoming (downwelling) shortwave solar radiation"),
    input(:surface_longwave_down, XY(), units=u"W/m^2", desc="Incoming (downwelling) longwave thermal radiation"),
    input(:daytime_length, XY(), default=12, units=u"hr", desc="Number of daytime hours varying with the season and orbital parameters")
)

"""
    shortwave_in(i, j, state, grid, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation})

Retrieve or compute the incoming/downwelling shortwave radiation at the current time step.
"""
shortwave_in(i, j, state, grid, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation}) where {PR} = state.surface_shortwave_down[i, j]

"""
    longwave_in(i, j, state, grid, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation})

Retrieve or compute the incoming/downwelling longwave radiation at the current time step.
"""
longwave_in(i, j, state, grid, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation}) where {PR} = state.surface_longwave_down[i, j]

"""
    daytime_length(i, j, state, grid, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation})

Retrieve the length of the day (in hours) at grid cell `i, j`. Defaults to a constant 12 hours if no input is provided.
"""
daytime_length(i, j, state, grid, ::AbstractAtmosphere{PR, <:LongShortWaveRadiation}) where {PR} = state.daytime_length[i, j]
