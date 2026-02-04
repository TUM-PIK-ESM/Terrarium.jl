@kwdef struct CoupledSoilEnergyModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Atmosphere<:AbstractAtmosphere,
    SurfaceEnergy<:AbstractSurfaceEnergyBalance,
    SoilProcesses<:AbstractSoil{NF},
    Initializer<:AbstractInitializer,
} <: AbstractSoilModel{NF, GridType}
    "Spatial discretization"
    grid::GridType

    "Near-surface atmospheric conditions"
    atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Surface energy balance"
    surface_energy_balance::SurfaceEnergy = SurfaceEnergyBalance(eltype(grid))

    "Soil processes"
    soil::SoilProcesses = SoilEnergyWaterCarbon(eltype(grid))

    "Initializer for coupled model"
    initializer::Initializer = DefaultInitializer()
end

function processes(model::CoupledSoilEnergyModel)
    return (
        model.atmosphere,
        model.surface_energy_balance,
        processes(model.soil)...,
    )
end

function variables(model::CoupledSoilEnergyModel)
    atmos_vars = variables(model.atmosphere)
    seb_vars = variables(model.surface_energy_balance)
    soil_vars = variables(model.soil)
    return (
        atmos_vars...,
        seb_vars...,
        soil_vars...,
    )
end

function initialize(
    model::CoupledSoilEnergyModel{NF};
    clock = Clock(time=zero(NF)),
    boundary_conditions = (;),
    fields = (;),
    input_variables = ()
) where {NF}
    grid = model.grid
    # Collect all model variables
    vars = Variables(variables(model)..., input_variables...)
    # Create ground heat flux field for coupling
    ground_heat_flux = initialize(vars.auxiliary.ground_heat_flux, grid, clock)
    ground_heat_flux_bc = GroundHeatFlux(ground_heat_flux)
    # Merge with user-specified BCs and fields
    boundary_conditions = merge_boundary_conditions(boundary_conditions, ground_heat_flux_bc)
    fields = (; ground_heat_flux, fields...)
    # Return the initialized model state
    return initialize(vars, grid; clock, boundary_conditions, fields)
end

function initialize!(state, model::CoupledSoilEnergyModel)
    # Call initialize! with model initializer
    initialize!(state, model, model.initializer)
    # Then for soil model
    initialize!(state, model.grid, model.soil)
end

function compute_auxiliary!(state, model::CoupledSoilEnergyModel)
    grid = get_grid(model)
    atmosphere = model.atmosphere
    seb = model.surface_energy_balance
    constants = model.constants
    compute_auxiliary!(state, grid, atmosphere)
    compute_auxiliary!(state, grid, seb, atmosphere, constants)
    compute_auxiliary!(state, grid, soil, constants)
end

function compute_tendencies!(state, model::CoupledSoilEnergyModel)
    grid = get_grid(model)
    constants = model.constants
    # Compute soil tendencies
    compute_tendencies!(state, grid, soil, constants)
end
