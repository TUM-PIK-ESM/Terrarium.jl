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