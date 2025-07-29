# land surface model types
abstract type AbstractSnowModel <: AbstractLandModel end
abstract type AbstractVegetationModel <: AbstractLandModel end
abstract type AbstractEnergyBalanceModel <: AbstractLandModel end
abstract type AbstractHydrologyModel <: AbstractLandModel end

struct LandSurfaceModel{
    Grid<:AbstractLandGrid,
    SnowModel<:AbstractSnowModel,
    VegetationModel<:AbstractVegetationModel,
    EnergyBalanceModel<:AbstractEnergyBalanceModel,
    HydrologyModel<:AbstractHydrologyModel,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractLandSurfaceModel
    "Spatial grid"
    grid::Grid

    "Snow scheme"
    snow::SnowModel

    "Vegetation dynamics"
    vegetation::VegetationModel

    "Surface energy balance"
    energy::EnergyBalanceModel
    
    "Hydrology type"
    hydrology::HydrologyModel

    "Time stepping scheme"
    time_stepping::TimeStepper
end
