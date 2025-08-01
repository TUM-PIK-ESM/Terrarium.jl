@kwdef struct SoilModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Stratigraphy<:AbstractStratigraphy,
    SoilEnergy<:AbstractSoilEnergyBalance,
    SoilHydrology<:AbstractSoilHydrology,
    Biogeochemistry<:AbstractSoilBiogeochemistry,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractSoilModel
    "Spatial grid type"
    grid::GridType

    "Stratigraphy of the soil"
    strat::Stratigraphy = HomogeneousSoil()

    "Soil energy balance"
    energy::SoilEnergy = SoilEnergyBalance()

    "Soil hydrology/water balance"
    hydrology::SoilHydrology = ImmobileSoilWater()

    "Soil biogeochemistry"
    biogeochem::Biogeochemistry = ConstantSoilCarbonDenisty()

    "Physical constants"
    constants::PhysicalConstants{NF} = PhysicalConstants()

    "Boundary conditions"
    boundary_conditions::BoundaryConditions = PrescribedFluxes()

    "State variable initializer"
    initializer::Initializer = FieldInitializer()

    "Timestepping type"
    time_stepping::TimeStepper = ForwardEuler()
end

# SoilModel getter methods

get_stratigraphy(model::SoilModel) = model.strat

get_energy_balance(model::SoilModel) = model.energy

get_hydrology(model::SoilModel) = model.hydrology

get_biogeochemistry(model::SoilModel) = model.biogeochem

get_constants(model::SoilModel) = model.constants

# Forwarded methods

freezecurve(model::SoilModel) = freezecurve(model.strat)

# Model interface methods

function variables(model::SoilModel)
    strat_vars = variables(model.strat)
    hydrology_vars = variables(model.energy)
    energy_vars = variables(model.energy)
    bgc_vars = variables(model.biogeochem)
    # combine all variables into one tuple
    return (strat_vars..., hydrology_vars..., energy_vars..., bgc_vars...)
end

function initialize!(state, model::SoilModel)
    # TODO
end

function update_state!(state, model::SoilModel)
    grid = get_grid(model)
    launch!(grid, _update_state_soil_kernel, state, model)
    return nothing
end

function compute_tendencies!(state, model::SoilModel)
    fill_halo_regions!(state.temperature)
    grid = get_grid(model)
    launch!(grid, _compute_tendencies_soil_kernel, state, model)
    return nothing
end

function timestep!(state, model::SoilModel, euler::ForwardEuler, dt=get_dt(euler))
    update_state!(state, model)
    compute_tendencies!(state, model)
    # just timestep energy state for now;
    # ideally, the timestepper should do this automatically for all variables
    @. state.enthalpy += dt*state.enthalpy_tendency
    return nothing
end

@kernel function _update_state_soil_kernel(state, model::SoilModel)
    i, j, k = @index(Global, NTuple)
    update_state!(i, j, k, state, model, model.strat)
    update_state!(i, j, k, state, model, model.hydrology)
    update_state!(i, j, k, state, model, model.energy)
    update_state!(i, j, k, state, model, model.biogeochem)
end

@kernel function _compute_tendencies_soil_kernel(state, model::SoilModel)
    i, j, k = @index(Global, NTuple)
    compute_tendencies!(i, j, k, state, model, model.strat)
    compute_tendencies!(i, j, k, state, model, model.hydrology)
    compute_tendencies!(i, j, k, state, model, model.energy)
    compute_tendencies!(i, j, k, state, model, model.biogeochem)
end

@inline function soil_characteristic_fractions(i, j, k, state, model)
    sat = state.pore_water_ice_saturation[i, j, k]
    por = porosity(i, j, k, state, model, get_stratigraphy(model))
    ## there is some slight redundant computation here; consider merging into one method?
    org = organic_fraction(i, j, k, state, model, get_biogeochemistry(model))
    return (; sat, por, org)
end

@inline function soil_volumetric_fractions(i, j, k, state, model)
    # get characteristic fractions
    sat, por, org = soil_characteristic_fractions(i, j, k, state, model)
    # get fraction of unfrozen pore water
    liq = state.liquid_water_fraction[i, j, k]
    # calculate volumetric fractions
    water_ice = sat*por
    water = water_ice*liq
    ice = water_ice*(1-liq)
    air = (1-sat)*por
    mineral = (1-por)*(1-org)
    organic = (1-por)*org
    return (; water, ice, air, mineral, organic)
end
