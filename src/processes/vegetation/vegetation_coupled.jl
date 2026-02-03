@kwdef struct DynamicVegetation{
    NF,
    Photosynthesis <: AbstractPhotosynthesis,
    StomatalConducatance <: AbstractStomatalConductance,
    AutotrophicRespiration <: AbstractAutotrophicRespiration,
    Phenology <: AbstractPhenology,
    CarbonDynamics <: AbstractVegetationCarbonDynamics,
    VegetationDynamics <: AbstractVegetationDynamics,
    RootDistribution <: AbstractRootDistribution
} <: AbstractVegetation{NF}
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
end

function compute_auxiliary!(
    state, grid,
    veg::DynamicVegetation,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants,
    args...
)
    # Compute auxiliary variables for each component
    # Veg. carbon dynamics: needs C_veg(t-1) and computes LAI_b(t-1)
    compute_auxiliary!(state, grid, veg.carbon_dynamics)

    # Phenology: needs LAI_b(t-1) and computes LAI(t-1) and phen(t-1)
    compute_auxiliary!(state, grid, veg.phenology)

    # Stomatal conductance: needs atm. inputs(t) and computes λc(t)
    # TODO: Note the (implicit) circular dependency between photosynthesis and stomatal conductance;
    # can this be refactored?
    compute_auxiliary!(state, grid, veg.stomatal_conductance, veg.photosynthesis, atmos, constants)

    # Photosynthesis: needs atm. inputs(t), λc(t), LAI(t-1), and computes Rd(t) and GPP(t)
    compute_auxiliary!(state, grid, veg.photosynthesis, veg.stomatal_conductance, atmos)

    # Autotrophic respiration: needs atm. inputs(t), GPP(t), Rd(t), C_veg(t-1), phen(t-1) and computes Ra(t) and NPP(t)
    compute_auxiliary!(state, grid, veg.autotrophic_respiration, veg.carbon_dynamics, atmos)
    
    # Note: vegetation_dynamics compute_auxiliary! does nothing for now
    compute_auxiliary!(state, grid, veg.vegetation_dynamics)
    return nothing
end

function compute_tendencies!(state, veg::DynamicVegetation, args...)
    # Needs NPP(t), C_veg(t-1), LAI_b(t-1) and computes tendency for C_veg
    compute_tendencies!(state, grid, veg.carbon_dynamics)

    # Needs NPP(t), C_veg(t-1), LAI_b(t-1), ν(t-1) and computes tendency for ν
    compute_tendencies!(state, grid, veg.vegetation_dynamics, veg.carbon_dynamics)

    return nothing
end
