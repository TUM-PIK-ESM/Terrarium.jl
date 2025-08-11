# Initial concept of what a semi-complete land model might look like.
struct LandModel{
    NF,
    GridType<:AbstractLandGrid,
    GroundModel<:AbstractGroundModel,
    SnowModel<:AbstractSnowModel,
    VegetationModel<:AbstractVegetationModel,
    EnergyBalanceModel<:AbstractEnergyBalanceModel,
    HydrologyModel<:AbstractHydrologyModel,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractLandModel{NF}
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
    
    "Surface hydrology model"
    hydrology::HydrologyModel
    
    "Time stepping scheme"
    time_stepping::TimeStepper
end
