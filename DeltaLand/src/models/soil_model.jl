struct SoilModel{
    GridType<:AbstractLandGrid,
    Stratigraphy<:AbstractStratigraphy,
    SoilEnergyBalance<:AbstractSoilEnergyBalance,
    SoilHydrology<:AbstractSoilHydrology,
    Biogeochemistry<:AbstractSoilBiogeochemistry,
    TimeStepper<:AbstractTimeStepper,
} <: AbstractSoilModel
    "Spatial grid type"
    grid::GridType

    "Stratigraphy of the soil"
    strat::Stratigraphy

    "Soil energy balance"
    energy::SoilEnergyBalance

    "Soil hydrology/water balance"
    hydrology::SoilHydrology

    "Soil biogeochemistry"
    biogeochem::Biogeochemistry

    "Timestepping type"
    timestepping::TimeStepper
end

function initialize(model::SoilModel, initializer; sim_kwargs...)
    # Extract grid from model
    grid = get_grid(model)

    # Initialize state variables
    # TODO: this is a temporary solution until we find a better solution;
    # Ideally we would like these variables to be declared by the individual process
    # implementations and then merged somehow, similar to what CryoGrid.jl does.
    enthalpy = initialize3D(initializer, grid, :enthalpy)
    enthalpy_tendency = initialize3D(initializer, grid, :enthalpy_tendency)
    temperature = initialize3D(initializer, grid, :temperature)
    water_ice_saturation = initialize3D(initializer, grid, :water_ice_saturation)
    water_saturation = initialize3D(initializer, grid, :water_ice_saturation)
    # 2D boundary variables
    ground_heat_flux = initialize2D(initializer, grid, :ground_heat_flux)
    geothermal_heat_flux = initialize2D(initializer, grid, :geothermal_heat_flux)
    
    # Create model state and simulation
    # TODO: it would be cool if we could formally specify constitutive relations like temperature <--> enthalpy
    state = StateVariables(
        prognostic=(enthalpy,),
        tendencies=(enthalpy_tendency,),
        auxiliary=(temperature, water_ice_saturation, water_saturation, ground_heat_flux, geothermal_heat_flux)
    )
    return Simulation(model, state, initializer.clock)
end

function update_state!(state, model::SoilModel)
    fill_halo_regions!(state)
    grid = get_grid(model)
    # TODO: shouldn't the workspec be determined by the grid type?
    launch!(grid.architecture, grid, :xyz, _update_state_soil_kernel, state, model)
    return nothing
end

function compute_tendencies!(state, model::SoilModel)
    grid = get_grid(model)
    # TODO: shouldn't the workspec be determined by the grid type?
    launch!(grid.architecture, grid, :xyz, _compute_tendencies_soil_kernel, state, model)
    return nothing
end

function timestep!(state, model::SoilModel, euler::ForwardEuler, dt=get_dt(euler))
    update_state!(state, model)
    compute_tendencies!(state, model)
    # just timestep energy state for now;
    # ideally, the timestepper should do this automatically for all variables
    state.enthalpy .+= dt*state.enthalpy_tendency
    return nothing
end

@kernel function _update_state_soil_kernel(state, model::SoilModel)
    i, j, k = @index(Global, NTuple)
    update_state(i, j, k, state, model, model.strat)
    update_state(i, j, k, state, model, model.hydrology)
    update_state(i, j, k, state, model, model.energy)
    update_state(i, j, k, state, model, model.biogeochem)
end

@kernel function _compute_tendencies_soil_kernel(state, model::SoilModel)
    i, j, k = @index(Global, NTuple)
    compute_tendencies(i, j, k, state, model, model.strat)
    compute_tendencies(i, j, k, state, model, model.hydrology)
    compute_tendencies(i, j, k, state, model, model.energy)
    compute_tendencies(i, j, k, state, model, model.biogeochem)
end
