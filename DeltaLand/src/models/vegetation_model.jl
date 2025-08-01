"""
    $TYPEDEF

Model for natural (unmanaged) vegetation processes for a single plant functional type (PFT).
Multiple PFTs can be later handled with a `TiledVegetationModel` type that composes multiple
`VegetationModel`s with different parameters for each PFT.
"""
@kwdef struct VegetationModel{
    NF,
    Photosynthesis<:AbstractPhotosynthesis,
    StomatalConducatance<:AbstractStomatalConductance,
    AutotrophicRespiration<:AbstractAutotrophicRespiration,
    CarbonDynamics<:AbstractVegetationCarbonDynamics,
    VegetationDynamics<:AbstractVegetationDynamics,
    Phenology<:AbstractPhenology,
    GridType<:AbstractLandGrid{NF},
    BoundaryConditions<:AbstractBoundaryConditions,
    Initializer<:AbstractInitializer,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractSoilModel
    "Spatial grid type"
    grid::GridType

    "Photosynthesis scheme"
    photosynthesis::Photosynthesis # not prognostic

    "Stomatal conducantance scheme"
    stomatal_conductance::StomatalConducatance # not prognostic

    "Autotrophic respiration scheme"
    autotrophic_respiration::AutotrophicRespiration # not prognostic

    "Phenology scheme"
    phenology::Phenology # not prognostic

    "Vegetation carbon pool dynamics"
    carbon_dynamics::CarbonDynamics # prognostic

    "Vegetation population density or coverage fraction dynamics"
    vegetation_dynamics::VegetationDynamics # prognostic

    "Boundary conditions"
    boundary_conditions::BoundaryConditions = FieldBoundaryConditions()

    "State variable initializer"
    initializer::Initializer = FieldInitializer()

    "Timestepping type"
    time_stepping::TimeStepper = ForwardEuler()
end

variables(model::VegetationModel) = (
    variables(model.photosynthesis)...,
    variables(model.stomatal_conductance)...,
    variables(model.autotrophic_respiration)...,
    variables(model.phenology)...,
    variables(model.carbon_dynamics)...,
    variables(model.vegetation_dynamics)...,
)

function update_state!(state, model::VegetationModel)

end

function compute_tendencies!(state, model::VegetationModel)

end
