# Initial concept of what a semi-complete land model might look like.
struct LandModel{
    GridType<:AbstractLandGrid,
    LandSurfaceModel<:AbstractLandSurfaceModel,
    GroundModel<:AbstractGroundModel,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractLandModel
    "Spatial grid"
    grid::GridType

    "Surface energy, hydrology, and vegetation model"
    surface::LandSurfaceModel

    "Ground model"
    ground::GroundModel
    
    "Time stepping scheme"
    time_stepping::TimeStepper
end
