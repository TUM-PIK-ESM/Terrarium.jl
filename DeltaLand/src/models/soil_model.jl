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
} <: AbstractSoilModel
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
    boundary_conditions::BoundaryConditions = PrescribedFluxes()

    "State variable initializer"
    initializer::Initializer = FieldInitializers()

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
    strat_vars = variables(model.strat)
    hydrology_vars = variables(model.hydrology)
    energy_vars = variables(model.energy)
    bgc_vars = variables(model.biogeochem)
    # combine all variables into one tuple
    return tuplejoin(strat_vars, hydrology_vars, energy_vars, bgc_vars)
end

function compute_auxiliary!(state, model::SoilModel)
    grid = get_grid(model)
    launch!(grid, _compute_auxiliary_soil_kernel, state, model)
    return nothing
end

function compute_tendencies!(state, model::SoilModel)
    fill_halo_regions!(state.temperature)
    grid = get_grid(model)
    launch!(grid, _compute_tendencies_soil_kernel, state, model)
    return nothing
end

function timestep!(state, model::SoilModel, dt=get_dt(timestepper))
    compute_auxiliary!(state, model)
    compute_tendencies!(state, model)
    grid = get_grid(model)
    launch!(grid, _timestep_kernel, state, model, get_time_stepping(model), dt)
    return nothing
end

@kernel function _compute_auxiliary_soil_kernel(state, model::SoilModel)
    idx = @index(Global, NTuple)
    compute_auxiliary!(idx, state, model, model.strat)
    compute_auxiliary!(idx, state, model, model.hydrology)
    compute_auxiliary!(idx, state, model, model.energy)
    compute_auxiliary!(idx, state, model, model.biogeochem)
end

@kernel function _compute_tendencies_soil_kernel(state, model::SoilModel)
    idx = @index(Global, NTuple)
    compute_tendencies!(idx, state, model, model.strat)
    compute_tendencies!(idx, state, model, model.hydrology)
    compute_tendencies!(idx, state, model, model.energy)
    compute_tendencies!(idx, state, model, model.biogeochem)
end

@kernel function _timestep_kernel(state, model::SoilModel, euler::ForwardEuler, dt)
    idx = @index(Global, NTuple)
    i, j, k = idx
    # timestep for internal energy
    state.internal_energy[i, j, k] = state.internal_energy[i, j, k] + dt*state.internal_energy_tendency[i, j, k]
    # apply inverse closure relation to update temperature
    energy_to_temperature!(idx, state, model, get_freezecurve(get_soil_hydrology(model)))
end

@inline function soil_characteristic_fractions(idx, state, model::SoilModel)
    sat = state.pore_water_ice_saturation[idx...]
    por = porosity(idx, state, model, get_soil_hydrology(model))
    ## there is some slight redundant computation here; consider merging into one method?
    org = organic_fraction(idx, state, model, get_biogeochemistry(model))
    return (; sat, por, org)
end

@inline function soil_volumetric_fractions(idx, state, model::SoilModel)
    # get characteristic fractions
    (; sat, por, org) = soil_characteristic_fractions(idx, state, model)
    # get fraction of unfrozen pore water
    liq = state.liquid_water_fraction[idx...]
    # calculate volumetric fractions
    water_ice = sat*por
    water = water_ice*liq
    ice = water_ice*(1-liq)
    air = (1-sat)*por
    mineral = (1-por)*(1-org)
    organic = (1-por)*org
    return (; water, ice, air, mineral, organic)
end

# Initialization

@kernel function _initialize_kernel(state, model::SoilModel, initializer)
    idx = @index(Global, NTuple)
    initialize!(idx, state, model, initializer)
end

function initialize!(state, model::SoilModel, initializer::AbstractInitializer)
    grid = get_grid(model)
    launch!(grid, _initialize_kernel, state, model, initializer)
end

function initialize!(idx, state, model::SoilModel, initializer::AbstractInitializer)
    # TODO: need a more comprehensive initialization scheme for all soil model components
    hydrology = get_soil_hydrology(model)
    # Note that this assumes temperature has already been iniitialized
    temperature_to_energy!(idx, state, model, get_freezecurve(hydrology))
end
