"""
    $TYPEDEF

Model for natural (unmanaged) vegetation processes for a single plant functional type (PFT).
Multiple PFTs can be later handled with a `TiledVegetationModel` type that composes multiple
`VegetationModel`s with different parameters for each PFT.

Properties:
$TYPEDFIELDS
"""
@kwdef struct VegetationModel{
    NF,
    Photosynthesis<:AbstractPhotosynthesis,
    StomatalConducatance<:AbstractStomatalConductance,
    AutotrophicRespiration<:AbstractAutotrophicRespiration,
    Phenology<:AbstractPhenology,
    CarbonDynamics<:AbstractVegetationCarbonDynamics,
    VegetationDynamics<:AbstractVegetationDynamics,
    RootDistribution<:AbstractRootDistribution,
    Atmosphere<:AbstractAtmosphere,
    GridType<:AbstractLandGrid{NF},
    Constants<:PhysicalConstants{NF},
    Initializer<:AbstractInitializer,
} <: AbstractVegetationModel{NF, GridType}
    "Spatial grid type"
    grid::GridType

    "Atmospheric input configuration"
    atmosphere::Atmosphere

    "Photosynthesis scheme"
    photosynthesis::Photosynthesis = LUEPhotosynthesis(eltype(grid)) # not prognostic

    "Stomatal conducantance scheme"
    stomatal_conductance::StomatalConducatance = MedlynStomatalConductance(eltype(grid)) # not prognostic

    "Autotrophic respiration scheme"
    autotrophic_respiration::AutotrophicRespiration = PALADYNAutotrophicRespiration(eltype(grid)) # not prognostic

    "Phenology scheme"
    phenology::Phenology = PALADYNPhenology(eltype(grid)) # not prognostic

    "Vegetation carbon pool dynamics"
    carbon_dynamics::CarbonDynamics = PALADYNCarbonDynamics(eltype(grid)) # prognostic

    "Vegetation population density or coverage fraction dynamics"
    vegetation_dynamics::VegetationDynamics =  PALADYNVegetationDynamics(eltype(grid)) # prognostic

    "Vegetation root distribution"
    root_distribution::RootDistribution = StaticExponentialRootDistribution(eltype(grid))

    "Physical constants"
    constants::Constants = PhysicalConstants(eltype(grid))

    "State variable initializer"
    initializer::Initializer = DefaultInitializer()
end

# VegetationModel getter methods
get_atmosphere(model::VegetationModel) = model.atmosphere

get_photosynthesis(model::VegetationModel) = model.photosynthesis

get_stomatal_conductance(model::VegetationModel) = model.stomatal_conductance

get_autotrophic_respiration(model::VegetationModel) = model.autotrophic_respiration

get_phenology(model::VegetationModel) = model.phenology

get_vegetation_carbon_dynamics(model::VegetationModel) = model.carbon_dynamics

get_vegetation_dynamics(model::VegetationModel) = model.vegetation_dynamics

get_root_distribution(model::VegetationModel) = model.root_distribution

get_constants(model::VegetationModel) = model.constants

# Model interface methods
variables(model::VegetationModel) = tuplejoin(
    variables(model.atmosphere),
    variables(model.photosynthesis),
    variables(model.stomatal_conductance),
    variables(model.autotrophic_respiration),
    variables(model.phenology),
    variables(model.carbon_dynamics),
    variables(model.vegetation_dynamics),
    variables(model.root_distribution)
)

processes(model::VegetationModel) = (
    model.atmosphere,
    model.photosynthesis,
    model.stomatal_conductance,
    model.autotrophic_respiration,
    model.carbon_dynamics,
    model.vegetation_dynamics
)

function compute_auxiliary!(state, model::VegetationModel)
    # Compute auxiliary variables for each component
    # Veg. carbon dynamics: needs C_veg(t-1) and computes LAI_b(t-1)
    compute_auxiliary!(state, model, model.carbon_dynamics)

    # Phenology: needs LAI_b(t-1) and computes LAI(t-1) and phen(t-1)
    compute_auxiliary!(state, model, model.phenology)

    # Stomatal conductance: needs atm. inputs(t) and computes λc(t)
    compute_auxiliary!(state, model, model.stomatal_conductance)

    # Photosynthesis: needs atm. inputs(t), λc(t), LAI(t-1), and computes Rd(t) and GPP(t)
    compute_auxiliary!(state, model, model.photosynthesis)

    # Autotrophic respiration: needs atm. inputs(t), GPP(t), Rd(t), C_veg(t-1), phen(t-1) and computes Ra(t) and NPP(t)
    compute_auxiliary!(state, model, model.autotrophic_respiration)
    
    # Note: vegetation_dynamics compute_auxiliary! does nothing for now
    return nothing
end

function compute_tendencies!(state, model::VegetationModel)
    # Needs NPP(t), C_veg(t-1), LAI_b(t-1) and computes tendency for C_veg
    compute_tendencies!(state, model, model.carbon_dynamics)

    # Needs NPP(t), C_veg(t-1), LAI_b(t-1), ν(t-1) and computes tendency for ν
    compute_tendencies!(state, model, model.vegetation_dynamics)

    return nothing
end
