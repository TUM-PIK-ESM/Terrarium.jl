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
    CarbonDynamics<:AbstractVegetationCarbonDynamics,
    VegetationDynamics<:AbstractVegetationDynamics,
    Phenology<:AbstractPhenology,
    GridType<:AbstractLandGrid{NF},
    Constants<:PhysicalConstants{NF},
    BoundaryConditions<:AbstractBoundaryConditions,
    Initializer<:AbstractInitializer,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractVegetationModel{NF}
    "Spatial grid type"
    grid::GridType

    "Photosynthesis scheme"
    photosynthesis::Photosynthesis = LUEPhotosynthesis() # not prognostic

    "Stomatal conducantance scheme"
    stomatal_conductance::StomatalConducatance = MedlynStomatalConductance() # not prognostic

    "Autotrophic respiration scheme"
    autotrophic_respiration::AutotrophicRespiration = PALADYNAutotrophicRespiration() # not prognostic

    "Phenology scheme"
    phenology::Phenology = PALADYNPhenology() # not prognostic

    "Vegetation carbon pool dynamics"
    carbon_dynamics::CarbonDynamics = PALADYNCarbonDynamics() # prognostic

    "Vegetation population density or coverage fraction dynamics"
    vegetation_dynamics::VegetationDynamics =  PALADYNVegetationDynamics() # prognostic

    "Physical constants"
    constants::Constants = PhysicalConstants{eltype(grid)}()

    "Boundary conditions"
    boundary_conditions::BoundaryConditions = DefaultBoundaryConditions()

    "State variable initializer"
    initializer::Initializer = Initializers()

    "Timestepping scheme"
    time_stepping::TimeStepper = ForwardEuler()
end

# VegetationModel getter methods
get_photosynthesis(model::VegetationModel) = model.photosynthesis

get_stomatal_conductance(model::VegetationModel) = model.stomatal_conductance

get_autotrophic_respiration(model::VegetationModel) = model.autotrophic_respiration

get_phenology(model::VegetationModel) = model.phenology

get_carbon_dynamics(model::VegetationModel) = model.carbon_dynamics

get_vegetation_dynamics(model::VegetationModel) = model.vegetation_dynamics

get_constants(model::VegetationModel) = model.constants


# Model interface methods
variables(model::VegetationModel) = (
    variables(model.photosynthesis)...,
    variables(model.stomatal_conductance)...,
    variables(model.autotrophic_respiration)...,
    variables(model.phenology)...,
    variables(model.carbon_dynamics)...,
    variables(model.vegetation_dynamics)...,
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
    # Fill halo regions for fields with boundary conditions
    fill_halo_regions!(state)

    # Needs NPP(t), C_veg(t-1), LAI_b(t-1) and computes C_veg_tendency
    compute_tendencies!(state, model, model.carbon_dynamics)

    # Needs NPP(t), C_veg(t-1), LAI_b(t-1), ν(t-1) and computes ν_tendency
    compute_tendencies!(state, model, model.vegetation_dynamics)

    return nothing
end

function timestep!(state, model::VegetationModel, euler::ForwardEuler, dt=get_dt(timestepper))
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    grid = get_grid(model)
    launch!(grid, :xy, timestep_vegetation_kernel!, state, euler, dt)
    return nothing
end

@kernel function timestep_vegetation_kernel!(state, euler::ForwardEuler, dt)
    i, j = @index(Global, NTuple)
    # Update vegetation carbon pool, compute C_veg(t)
    state.C_veg[i, j] = step(euler, state.C_veg[i, j], state.C_veg_tendency[i, j], dt)
    # Update vegetation fraction, compute ν(t)
    state.ν[i, j] = step(euler, state.ν[i, j], state.ν_tendency[i, j], dt)
end
