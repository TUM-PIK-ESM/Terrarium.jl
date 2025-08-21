@kwdef struct SoilInitializer{
    EnergyInit,
    HydrologyInit,
    StratInit,
    BGCInit
} <: AbstractInitializer
    "Soil energy/temperature state initializer"
    energy::EnergyInit = DefaultInitializer()

    "Soil hydrology state initializer"
    hydrology::HydrologyInit = DefaultInitializer()

    "Soil stratigraphy state initializer"
    strat::StratInit = DefaultInitializer()

    "Soil biogeochemistry state initializer"
    biogeochem::BGCInit = DefaultInitializer()
end

function initialize!(state, model::AbstractModel, init::SoilInitializer)
    initialize!(state, model, init.energy)
    initialize!(state, model, init.hydrology)
    initialize!(state, model, init.strat)
    initialize!(state, model, init.biogeochem)
end

function get_field_initializers(inits::SoilInitializer)
    return merge(
        get_field_initializers(inits.energy),
        get_field_initializers(inits.hydrology),
        get_field_initializers(inits.strat),
        get_field_initializers(inits.biogeochem),
    )
end

# Soil energy initializers

"""
Creates a constant soil temperature initializer.
"""
ConstantInitialSoilTemperature(T₀) = Initializers(temperature = T₀)

# TODO: Add "real" thermal steady state Initializer
"""
Computes a linear temperature profile in quasi-steady state based on the given
surface temperature, geothermal heat flux, and bulk thermal conductivity.
"""
QuasiThermalSteadyState(T₀, Qgeo, k_eff) = Initializers(temperature = (x,z) -> T₀ + Qgeo / k_eff*z)

"""
Creates a piecwise linear temperature initializer from the given knots.

```julia
initializer = PiecewiseLinearInitialSoilTemperature(
    0.0u"m" => 5.0, # always in °C!
    0.5u"m" => 2.0,
    1.0u"m" => 1.0,
    10.0u"m" => 1.5,
    ...
)
```
"""
function PiecewiseLinearInitialSoilTemperature(knots::Pair{<:LengthQuantity}...)
    f = piecewise_linear(knots...)
    return Initializers(temperature = (x, z) -> f(z))
end
