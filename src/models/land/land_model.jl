# Initial concept of what a semi-complete land model might look like.
struct LandModel{
    NF,
    GridType<:AbstractLandGrid,
    GroundModel<:AbstractGroundModel,
    SnowModel<:AbstractSnowModel,
    VegetationModel<:AbstractVegetationModel,
    SurfaceEnergyBalanceModel<:AbstractSurfaceEnergyBalanceModel,
    HydrologyModel<:AbstractHydrologyModel,
    BoundaryConditions<:AbstractBoundaryConditions,
    Initializer<:AbstractInitializer,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractLandModel{NF, GridType, TimeStepper}
    "Spatial grid"
    grid::GridType

    "Ground model"
    ground::GroundModel

    "Snow scheme"
    snow::SnowModel

    "Vegetation dynamics"
    vegetation::VegetationModel

    "Surface energy balance"
    energy::SurfaceEnergyBalanceModel
    
    "Surface hydrology model"
    hydrology::HydrologyModel

    "Bounday conditions"
    boundary_conditions::BoundaryConditions

    "State variable initializer"
    initializer::Initializer
    
    "Time stepping scheme"
    time_stepping::TimeStepper
end

# TODO
variables(::LandModel) = ()

function compute_auxiliary!(state, ::LandModel)
    # TODO
end

function compute_tendencies!(state, ::LandModel)
    # TODO
end
