"""
    $TYPEDEF

General implementation of a 1D column model of soil energy, water, and carbon transport.

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct SoilModel{
    NF,
    GridType <: AbstractLandGrid{NF},
    Soil <: AbstractSoil{NF},
    Constants <: PhysicalConstants{NF},
    Initializer <: AbstractInitializer,
} <: AbstractSoilModel{NF, GridType}
    "Spatial grid type"
    grid::GridType

    "Soil processes"
    soil::Soil = SoilEnergyHydrologyBGC(eltype(grid))

    "Physical constants"
    constants::Constants = PhysicalConstants(eltype(grid))

    "State variable initializer"
    initializer::Initializer = SoilInitializer()
end

# Model interface methods

function initialize!(state, model::SoilModel)
    # run model/field initializers
    initialize!(state, model, model.initializer)
    # run process initializers
    initialize!(state, model.grid, model.soil, model.constants)
end

function compute_auxiliary!(state, model::SoilModel)
    compute_auxiliary!(state, model.grid, model.soil, model.constants)
    return nothing
end

function compute_tendencies!(state, model::SoilModel)
    compute_tendencies!(state, model.grid, model.soil, model.constants)
    return nothing
end

# Closures

function closure!(state, model::SoilModel)
    closure!(state, model.grid, model.soil, model.constants)
end

function invclosure!(state, model::SoilModel)
    invclosure!(state, model.grid, model.soil, model.constants)
end
