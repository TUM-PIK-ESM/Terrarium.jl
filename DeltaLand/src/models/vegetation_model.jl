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
} <: AbstractVegetationModel
    "Spatial grid type"
    grid::GridType

    "Photosynthesis scheme"
    photosynthesis::Photosynthesis = LUEPhotosynthesis() # not prognostic

    "Stomatal conducantance scheme"
    stomatal_conductance::StomatalConducatance = MedlynStomatalConductance() # not prognostic

    "Autotrophic respiration scheme"
    autotrophic_respiration::AutotrophicRespiration = PaladynAutotrophicRespiration() # not prognostic

    "Phenology scheme"
    phenology::Phenology = PaladynPhenology() # not prognostic

    "Vegetation carbon pool dynamics"
    carbon_dynamics::CarbonDynamics = PaladynCarbonDynamics() # prognostic

    "Vegetation population density or coverage fraction dynamics"
    vegetation_dynamics::VegetationDynamics =  PaladynVegetationDynamics() # prognostic

    "Physical constants"
    constants::Constants = PhysicalConstants{eltype(grid)}()

    "Boundary conditions"
    boundary_conditions::BoundaryConditions = FieldBoundaryConditions()

    "State variable initializer"
    initializer::Initializer = FieldInitializers()

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
    grid = get_grid(model)
    launch!(grid, _compute_auxiliary_vegetation_kernel, state, model)
    return nothing
end

function compute_tendencies!(state, model::VegetationModel)
    grid = get_grid(model)
    launch!(grid, _compute_tendencies_vegetation_kernel, state, model)
    return nothing
end

function timestep!(state, model::VegetationModel, dt=get_dt(timestepper))
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    grid = get_grid(model)
    launch!(grid, _timestep_vegetation_kernel, state, model, get_time_stepping(model), dt)
    return nothing
end

# Kernel functions
@kernel function _compute_auxiliary_vegetation_kernel(state, model::VegetationModel)
    idx = @index(Global, NTuple)

    # Compute auxiliary variables for each component
    compute_auxiliary!(idx, state, model, model.stomatal_conductance) # This computes λc
    compute_auxiliary!(idx, state, model, model.photosynthesis) # This computes GPP
    compute_auxiliary!(idx, state, model, model.autotrophic_respiration) # This computes Ra and NPP
    compute_auxiliary!(idx, state, model, model.carbon_dynamics) # This computes LAI_b 
    # compute_auxiliary!(idx, state, model, model.vegetation_dynamics) 
    compute_auxiliary!(idx, state, model, model.phenology) # This computes LAI

end

@kernel function _compute_tendencies_vegetation_kernel(state, model::VegetationModel)
    idx = @index(Global, NTuple)
    compute_tendencies!(idx, state, model, model.carbon_dynamics)
    compute_tendencies!(idx, state, model, model.vegetation_dynamics)
end

@kernel function _timestep_vegetation_kernel(state, model::VegetationModel, euler::ForwardEuler, dt)
    idx = @index(Global, NTuple)
    i, j = idx
    # Update vegetation carbon pool
    state.C_veg[i, j] = state.C_veg[i, j] + 
                                  dt * state.C_veg_tendency[i, j]
    # Update vegetation fraction
    state.ν[i, j] = state.ν[i, j] + 
                               dt * state.ν_tendency[i, j]
end

# Initialization
@kernel function _initialize_vegetation_kernel(state, model::VegetationModel, initializer)
    idx = @index(Global, NTuple)
    initialize!(idx, state, model, initializer)
end

function initialize!(state, model::VegetationModel, initializer::AbstractInitializer)
    grid = get_grid(model)
    launch!(grid, _initialize_vegetation_kernel, state, model, initializer)
end
