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

# TODO
variables(::LandModel) = ()

function compute_auxiliary!(state, ::LandModel)
    # TODO
end

function compute_tendencies!(state, ::LandModel)
    # TODO
end
