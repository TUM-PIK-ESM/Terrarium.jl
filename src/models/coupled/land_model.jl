# Initial concept of what a semi-complete land model might look like.
struct LandModel{
    NF,
    GridType<:AbstractLandGrid,
    Atmosphere<:AbstractAtmosphere,
    SEB<:SurfaceEnergyBalance,
    GroundModel<:AbstractGroundModel,
    SnowModel<:AbstractSnowModel,
    VegetationModel<:AbstractVegetationModel,
    HydrologyModel<:AbstractHydrologyModel,
    Initializer<:AbstractInitializer,
} <: AbstractLandModel{NF, GridType}
    "Spatial grid"
    grid::GridType

    "Atmospheric inputs"
    atmosphere::Atmosphere

    "Surface energy balance"
    suface_energy_balance::SurfaceEnergyBalance

    "Ground model"
    ground::GroundModel

    "Snow scheme"
    snow::SnowModel

    "Vegetation dynamics"
    vegetation::VegetationModel
    
    "Surface hydrology model"
    hydrology::HydrologyModel

    "State variable initializer"
    initializer::Initializer
end

variables(::LandModel) = ()

processes(model::LandModel) = (
    model.atmosphere,
    model.surface_energy_balance,
    processes(model.ground)...,
    processes(model.snow)...,
    processes(model.hydrology)...,
)

function compute_auxiliary!(state, ::LandModel)
    # TODO
end

function compute_tendencies!(state, ::LandModel)
    # TODO
end
