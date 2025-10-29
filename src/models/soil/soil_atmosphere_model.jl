struct CoupledSoilAtmosphereModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Atmosphere<:AbstractAtmosphere,
    SurfaceEnergy<:AbstractSurfaceEnergyBalance,
    TimeStepper<:AbstractTimeStepper,
    SoilModel<:AbstractSoilModel{NF, GridType, TimeStepper},
    Initializer<:AbstractInitializer,
} <: AbstractSoilModel{NF, GridType, TimeStepper}
    "Spatial discretization"
    grid::GridType

    "Near-surface atmospheric conditions"
    atmosphere::Atmosphere

    "Surface energy balance"
    surface_energy_balance::SurfaceEnergy

    "Soil model"
    soil::SoilModel

    "Physical constants"
    constants::PhysicalConstants{NF}

    "Initializer for coupled model"
    initializer::Initializer
end

function CoupledSoilAtmosphereModel(
    soil::SoilModel{NF};
    atmosphere::AbstractAtmosphere = PrescribedAtmosphere(eltype(soil.grid)),
    surface_energy_balance::SurfaceEnergyBalance = SurfaceEnergyBalance(eltype(soil.grid)),
    constants::PhysicalConstants = PhysicalConstants(NF),
    initializer::AbstractInitializer = DefaultInitializer()
) where {NF}
    surface_fluxes = SoilBC(energy=GroundHeatFlux(), hydrology=DefaultBoundaryConditions())
    coupled_soil_bcs = SoilBoundaryConditions(
        NF,
        top = surface_fluxes,
        bottom = DefaultBoundaryConditions() # temporary hack
    )
    soil = setproperties(soil, boundary_conditions=coupled_soil_bcs)
    return CoupledSoilAtmosphereModel(soil.grid, atmosphere, surface_energy_balance, soil, constants, initializer)
end

get_boundary_conditions(model::CoupledSoilAtmosphereModel) = get_boundary_conditions(model.soil)

get_time_stepping(model::CoupledSoilAtmosphereModel) = model.soil.time_stepping

get_soil_energy_balance(model) = model.soil.energy

get_soil_hydrology(model) = model.soil.hydrology

get_stratigraphy(model) = model.soil.strat

get_biogeochemistry(model) = model.soil.biogeochem

get_constants(model) = model.soil.constants

function variables(model::CoupledSoilAtmosphereModel)
    atmos_vars = variables(model.atmosphere)
    seb_vars = variables(model.surface_energy_balance)
    soil_vars = variables(model.soil)
    return (
        atmos_vars...,
        seb_vars...,
        soil_vars...,
        # redefine ground temperature as auxiliary variable
        auxiliary(:ground_temperature, XY(), units=u"°C", desc="Temperature of the uppermost ground or soil grid cell in °C")
    )
end

function initialize!(state, model::CoupledSoilAtmosphereModel)
    initialize!(state, model, model.initializer)
    initialize!(state, model.soil)
end

function compute_auxiliary!(state, model::CoupledSoilAtmosphereModel)
    grid = get_grid(model)
    # set ground temperature field to uppermost soil temperature
    # TODO: Ideally the ground temperature field should *just be* the view of the uppermost
    # soil temeprature, but this currently isn't possible.
    set!(state.ground_temperature, view(state.temperature, :, :, grid.grid.Nz))
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.surface_energy_balance)
    compute_auxiliary!(state, model.soil)
end

function compute_tendencies!(state, model::CoupledSoilAtmosphereModel)
    compute_tendencies!(state, model, model.atmosphere)
    compute_tendencies!(state, model, model.surface_energy_balance)
    compute_tendencies!(state, model.soil)
end

function timestep!(state, model::CoupledSoilAtmosphereModel)
    # Just forward call to soil model
    timestep!(state, model.soil)
end
