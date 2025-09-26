"""
    $TYPEDEF

General implementation of a 1D column model of soil energy, water, and carbon transport.

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct SoilModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Stratigraphy<:AbstractStratigraphy,
    SoilEnergy<:AbstractSoilEnergyBalance,
    SoilHydrology<:AbstractSoilHydrology,
    Biogeochemistry<:AbstractSoilBiogeochemistry,
    Constants<:PhysicalConstants{NF},
    BoundaryConditions<:AbstractBoundaryConditions,
    Initializer<:AbstractInitializer,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractSoilModel{NF}
    "Spatial grid type"
    grid::GridType

    "Stratigraphy of the soil"
    strat::Stratigraphy = HomogeneousSoil()

    "Soil energy balance"
    energy::SoilEnergy = SoilEnergyBalance()

    "Soil hydrology/water balance"
    hydrology::SoilHydrology = SoilHydrology()

    "Soil biogeochemistry"
    biogeochem::Biogeochemistry = ConstantSoilCarbonDenisty()

    "Physical constants"
    constants::Constants = PhysicalConstants{eltype(grid)}()

    "Boundary conditions"
    boundary_conditions::BoundaryConditions = SoilBoundaryConditions(grid)

    "State variable initializer"
    initializer::Initializer = SoilInitializer()

    "Timestepping scheme"
    time_stepping::TimeStepper = ForwardEuler()
end

# SoilModel getter methods

get_stratigraphy(model::SoilModel) = model.strat

get_soil_energy_balance(model::SoilModel) = model.energy

get_soil_hydrology(model::SoilModel) = model.hydrology

get_biogeochemistry(model::SoilModel) = model.biogeochem

get_constants(model::SoilModel) = model.constants

# Model interface methods

function variables(model::SoilModel)
    bc_vars = variables(model.boundary_conditions)
    strat_vars = variables(model.strat)
    hydrology_vars = variables(model.hydrology)
    energy_vars = variables(model.energy)
    bgc_vars = variables(model.biogeochem)
    # combine all variables into one tuple
    return tuplejoin(bc_vars, strat_vars, hydrology_vars, energy_vars, bgc_vars)
end

function compute_auxiliary!(state, model::SoilModel)
    compute_auxiliary!(state, model, model.boundary_conditions)
    compute_auxiliary!(state, model, model.strat)
    compute_auxiliary!(state, model, model.hydrology)
    compute_auxiliary!(state, model, model.energy)
    compute_auxiliary!(state, model, model.biogeochem)
    return nothing
end

function compute_tendencies!(state, model::SoilModel)
    # Fill halo regions for fields with boundary conditions
    fill_halo_regions!(state)
    # Default implementation forwards the method dispatch to processes in the order:
    # Stratigraphy -> Hydrology -> Energy -> Biogeochemistry
    compute_tendencies!(state, model, model.strat)
    compute_tendencies!(state, model, model.hydrology)
    compute_tendencies!(state, model, model.energy)
    compute_tendencies!(state, model, model.biogeochem)
    return nothing
end

function timestep!(state, model::SoilModel, euler::ForwardEuler, dt=get_dt(timestepper))
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    launch!(
        model.grid,
        :xyz,
        timestep_kernel!,
        state,
        euler,
        model.energy,
        model.hydrology,
        model.strat,
        model.biogeochem,
        model.constants,
        dt
    )
    return nothing
end

@kernel function timestep_kernel!(
    state,
    euler::ForwardEuler,
    energy::AbstractSoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
    dt
)
    idx = @index(Global, NTuple)
    i, j, k = idx
    # timestep for internal energy
    state.internal_energy[i, j, k] = step(euler, state.internal_energy[i, j, k], state.internal_energy_tendency[i, j, k], dt)
    # apply inverse closure relation to update temperature
    fc = get_freezecurve(hydrology)
    energy_to_temperature!(idx, state, fc, energy, hydrology, strat, bgc, constants)
end

# Initialization

function initialize!(state, model::SoilModel)
    # run model initializer
    initialize!(state, model, model.initializer)
    # launch kernel for generic initailization routine
    grid = get_grid(model)
    launch!(
        grid,
        :xyz,
        initialize_kernel!,
        state,
        model.energy,
        model.hydrology,
        model.strat,
        model.biogeochem,
        model.constants
    )
end

@kernel function initialize_kernel!(
    state,
    energy::AbstractSoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
)
    idx = @index(Global, NTuple)
    # TODO: need a more comprehensive initialization scheme for all soil model components
    # Note that this assumes temperature has already been iniitialized!
    fc = get_freezecurve(hydrology)
    temperature_to_energy!(idx, state, fc, energy, hydrology, strat, bgc, constants)
end
