"""
    $TYPEDEF

General implementation of a 1D column model of soil energy, water, and carbon transport.

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct SoilModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Stratigraphy<:AbstractStratigraphy,
    SoilEnergy<:AbstractSoilEnergyBalance,
    SoilHydrology<:AbstractSoilHydrology,
    Biogeochemistry<:AbstractSoilBiogeochemistry,
    Constants<:PhysicalConstants{NF},
    Initializer<:AbstractInitializer,
} <: AbstractSoilModel{NF, GridType}
    "Spatial grid type"
    grid::GridType

    "Stratigraphy of the soil"
    strat::Stratigraphy = HomogeneousSoil(eltype(grid))

    "Soil energy balance"
    energy::SoilEnergy = SoilEnergyBalance(eltype(grid))

    "Soil hydrology/water balance"
    hydrology::SoilHydrology = SoilHydrology(eltype(grid))

    "Soil biogeochemistry"
    biogeochem::Biogeochemistry = ConstantSoilCarbonDensity(eltype(grid))

    "Physical constants"
    constants::Constants = PhysicalConstants(eltype(grid))

    "State variable initializer"
    initializer::Initializer = SoilInitializer()
end

# SoilModel getter methods

get_soil_stratigraphy(model::SoilModel) = model.strat

get_soil_energy_balance(model::SoilModel) = model.energy

get_soil_hydrology(model::SoilModel) = model.hydrology

get_soil_biogeochemistry(model::SoilModel) = model.biogeochem

get_constants(model::SoilModel) = model.constants

# Model interface methods

variables(model::SoilModel) = tuplejoin(
    variables(model.hydrology),
    variables(model.energy),
    variables(model.strat),
    variables(model.biogeochem),
)

processes(model::SoilModel) = (
    model.strat,
    model.energy,
    model.hydrology,
    model.biogeochem
)

get_closures(model::SoilModel) = (
    get_closure(model.hydrology),
    get_closure(model.energy)
)

function compute_auxiliary!(state, model::SoilModel)
    compute_auxiliary!(state, model, model.biogeochem)
    compute_auxiliary!(state, model, model.strat)
    compute_auxiliary!(state, model, model.hydrology)
    compute_auxiliary!(state, model, model.energy)
    return nothing
end

function compute_tendencies!(state, model::SoilModel)
    compute_tendencies!(state, model, model.biogeochem)
    compute_tendencies!(state, model, model.strat)
    compute_tendencies!(state, model, model.hydrology)
    compute_tendencies!(state, model, model.energy)
    return nothing
end

# Initialization

function initialize!(state, model::SoilModel)
    # run model/field initializers
    initialize!(state, model, model.initializer)
    # run process initializers
    initialize!(state, model, model.strat)
    initialize!(state, model, model.biogeochem)
    initialize!(state, model, model.hydrology)
    initialize!(state, model, model.energy)
end
