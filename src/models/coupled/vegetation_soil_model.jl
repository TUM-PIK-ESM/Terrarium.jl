@kwdef struct VegetationSoilModel{
        NF,
        GridType <: AbstractLandGrid{NF},
        Atmosphere <: AbstractAtmosphere,
        SEB <: AbstractSurfaceEnergyBalance,
        Hydrology <: AbstractSurfaceHydrology,
        Vegetation <: AbstractVegetation{NF},
        PAW <: AbstractPlantAvailableWater,
        Soil <: AbstractSoil{NF},
        Initializer <: AbstractInitializer,
    } <: AbstractLandModel{NF, GridType}
    "Spatial discretization"
    grid::GridType

    "Near-surface atmospheric conditions"
    atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Surface energy balance"
    surface_energy_balance::SEB = SurfaceEnergyBalance(eltype(grid))

    "Surface hydrology scheme"
    surface_hydrology::Hydrology = SurfaceHydrology(eltype(grid))

    "Vegetation model"
    vegetation::Vegetation = VegetationCarbon(eltype(grid))

    "Plant available water"
    plant_available_water::PAW = FieldCapacityLimitedPAW(eltype(grid))

    "Soil processes"
    soil::Soil = SoilEnergyWaterCarbon(eltype(grid))

    "Physical constants"
    constants::PhysicalConstants{NF} = PhysicalConstants(eltype(grid))

    "Initializer for coupled model"
    initializer::Initializer = DefaultInitializer()
end

function initialize(
        model::VegetationSoilModel{NF};
        clock = Clock(time = zero(NF)),
        boundary_conditions = (;),
        fields = (;),
        input_variables = ()
    ) where {NF}
    grid = get_grid(model)
    # Create ground heat flux field for coupling
    vars = Variables(variables(model)..., input_variables...)
    # Initialize BC fields for coupling
    ground_heat_flux = initialize(vars.auxiliary.ground_heat_flux, grid, clock, boundary_conditions, fields)
    infiltration = initialize(vars.auxiliary.infiltration, grid, clock, boundary_conditions, fields)
    ground_heat_flux_bc = GroundHeatFlux(ground_heat_flux)
    # Note that the hydrology module computes infiltration as positive so we need to negate it here
    # since fluxes are by convention positive upwards
    infiltration_bc = InfiltrationFlux(-infiltration)
    bcs = merge_boundary_conditions(boundary_conditions, ground_heat_flux_bc, infiltration_bc)
    # Merge user-defined fields with BC fields
    fields = merge((; ground_heat_flux, infiltration), fields)
    # Return the initialized model state
    return initialize(vars, grid; clock, boundary_conditions = bcs, fields)
end

function initialize!(state, model::VegetationSoilModel)
    # model initializer
    initialize!(state, model, model.initializer)
    # process initializers
    grid = get_grid(model)
    initialize!(state, grid, model.surface_energy_balance)
    initialize!(state, grid, model.surface_hydrology)
    # TODO: change when refactoring model/process types
    initialize!(state, grid, model.vegetation, model.atmosphere, model.constants)
    return initialize!(state, grid, model.soil, model.constants)
end

function compute_auxiliary!(state, model::VegetationSoilModel)
    grid = get_grid(model)
    # Compute auxiliary variables for atmosphere, if any
    compute_auxiliary!(state, grid, model.atmosphere)
    # Compute auxiliary variables for soil
    compute_auxiliary!(state, grid, model.soil, model.constants)
    # Compute auxiliary variables for surface hydrology
    compute_auxiliary!(state, grid, model.surface_hydrology, model.atmosphere, model.constants, model.soil)
    # Compute auxiliary variables for surface energy balance
    compute_auxiliary!(state, grid, model.surface_energy_balance, model.atmosphere, model.constants, model.surface_hydrology)
    compute_auxiliary!(state, grid, model.plant_available_water, model.soil)
    compute_auxiliary!(state, grid, model.vegetation, model.atmosphere, model.constants)
    # recompute surface energy fluxes
    return compute_surface_energy_fluxes!(state, grid, model.surface_energy_balance, model.atmosphere, model.constants, model.surface_hydrology)
end

function compute_tendencies!(state, model::VegetationSoilModel)
    grid = get_grid(model)
    compute_tendencies!(state, grid, model.surface_hydrology)
    compute_tendencies!(state, grid, model.soil, model.constants)
    return compute_tendencies!(state, grid, model.vegetation)
end

function closure!(state, model::VegetationSoilModel)
    grid = get_grid(model)
    return closure!(state, grid, model.soil, model.constants)
end

function invclosure!(state, model::VegetationSoilModel)
    grid = get_grid(model)
    return invclosure!(state, grid, model.soil, model.constants)
end
