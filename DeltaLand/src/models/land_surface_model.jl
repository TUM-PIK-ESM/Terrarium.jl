# land surface model types
abstract type AbstractSnowModel <: AbstractLandModel end
abstract type AbstractVegetationModel <: AbstractLandModel end
abstract type AbstractEnergyBalanceModel <: AbstractLandModel end
abstract type AbstractHydrologyModel <: AbstractLandModel end

struct LandSurfaceModel{
    SnowModel<:AbstractSnowModel,
    VegetationModel<:AbstractVegetationModel,
    EnergyBalanceModel<:AbstractEnergyBalanceModel,
    HydrologyModel<:AbstractHydrologyModel,
} <: AbstractLandSurfaceModel
    "Snow scheme"
    snow::SnowModel

    "Vegetation dynamics"
    vegetation::VegetationModel

    "Surface energy balance"
    energy::EnergyBalanceModel
    
    "Hydrology type"
    hydrology::HydrologyModel
end
