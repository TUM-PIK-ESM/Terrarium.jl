@kwdef struct CoupledSoilEnergyModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Atmosphere<:AbstractAtmosphere,
    SurfaceEnergy<:AbstractSurfaceEnergyBalance,
    SoilModel<:AbstractSoilModel{NF, GridType},
    Initializer<:AbstractInitializer,
} <: AbstractSoilModel{NF, GridType}
    "Spatial discretization"
    grid::GridType

    "Near-surface atmospheric conditions"
    atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Surface energy balance"
    surface_energy_balance::SurfaceEnergy = SurfaceEnergyBalance(eltype(grid))

    "Soil model"
    soil::SoilModel = SoilModel(eltype(grid))

    "Initializer for coupled model"
    initializer::Initializer = DefaultInitializer()
end

get_soil_energy_balance(model::CoupledSoilEnergyModel) = model.soil.energy

get_soil_hydrology(model::CoupledSoilEnergyModel) = model.soil.hydrology

get_soil_stratigraphy(model::CoupledSoilEnergyModel) = model.soil.strat

get_soil_biogeochemistry(model::CoupledSoilEnergyModel) = model.soil.biogeochem

get_constants(model::CoupledSoilEnergyModel) = model.soil.constants

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
    external_variables = ()
) where {NF}
    grid = get_grid(model)
    vars = Variables(variables(model)..., external_variables...)
    ground_heat_flux = initialize(vars.auxiliary.ground_heat_flux, grid, clock)
    ground_heat_flux_bc = GroundHeatFlux(ground_heat_flux)
    bcs = merge_recursive(boundary_conditions, ground_heat_flux_bc)
    return initialize(vars, model.grid, clock, boundary_conditions=bcs)
end

function initialize!(state, model::CoupledSoilEnergyModel)
    # Call initialize! with model initializer
    initialize!(state, model, model.initializer)
    # Then for soil model
    initialize!(state, model.soil)
end

function compute_auxiliary!(state, model::CoupledSoilEnergyModel)
    grid = get_grid(model)
    # set ground temperature field to uppermost soil temperature
    # TODO: Probably the ground temperature field should be initialized as a
    # view of the uppermost soil layer temperature to save memory
    set!(state.ground_temperature, view(state.temperature, :, :, grid.grid.Nz))
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.surface_energy_balance)
    compute_auxiliary!(state, model.soil)
end

function compute_tendencies!(state, model::CoupledSoilEnergyModel)
    compute_tendencies!(state, model, model.surface_energy_balance)
    compute_tendencies!(state, model.soil)
end
