# Initial concept of what a semi-complete land model might look like.
struct LandModel{
    GridType<:AbstractLandGrid,
    GroundModel<:AbstractGroundModel,
    SnowModel<:AbstractSnowModel,
    VegetationModel<:AbstractVegetationModel,
    EnergyBalanceModel<:AbstractEnergyBalanceModel,
    HydrologyModel<:AbstractHydrologyModel,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractLandModel
    "Spatial grid"
    grid::GridType

    "Ground model"
    ground::GroundModel

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
