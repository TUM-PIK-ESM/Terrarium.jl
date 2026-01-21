@kwdef struct VegetationSoilModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Atmosphere<:AbstractAtmosphere,
    SEB<:AbstractSurfaceEnergyBalance,
    Hydrology<:AbstractSurfaceHydrology,
    Vegetation<:AbstractVegetationModel{NF, GridType},
    SoilModel<:AbstractSoilModel{NF, GridType},
    Initializer<:AbstractInitializer,
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
    vegetation::Vegetation = VegetationModel(grid; atmosphere)

    "Soil model"
    soil::SoilModel = SoilModel(grid)

    "Initializer for coupled model"
    initializer::Initializer = DefaultInitializer()
end

get_surface_energy_balance(model::VegetationSoilModel) = model.surface_energy

get_surface_hydrology(model::VegetationSoilModel) = model.surface_hydrology

get_soil_energy_balance(model::VegetationSoilModel) = model.soil.energy

get_soil_hydrology(model::VegetationSoilModel) = model.soil.hydrology

get_soil_stratigraphy(model::VegetationSoilModel) = model.soil.strat

get_soil_biogeochemistry(model::VegetationSoilModel) = model.soil.biogeochem

# just ignore that these are duplicated in veg model for now...
get_constants(model::VegetationSoilModel) = model.soil.constants

processes(model::VegetationSoilModel) = (
    model.atmosphere,
    model.surface_energy_balance,
    model.surface_hydrology,
    processes(model.vegetation)...,
    processes(model.soil)...
)

function initialize(
    model::VegetationSoilModel{NF};
    clock = Clock(time=zero(NF)),
    boundary_conditions = (;),
    fields = (;),
    external_variables = ()
) where {NF}
    grid = get_grid(model)
    vars = Variables(variables(model)..., external_variables...)
    ground_heat_flux = initialize(vars.auxiliary.ground_heat_flux, grid, clock)
    infiltration = initialize(vars.auxiliary.infiltration, grid, clock)
    ground_heat_flux_bc = GroundHeatFlux(ground_heat_flux)
    infiltration_bc = InfiltrationFlux(infiltration)
    bcs = merge_recursive(boundary_conditions, ground_heat_flux_bc, infiltration_bc)
    fields = merge(fields, (; ground_heat_flux, infiltration))
    return initialize(vars, model.grid, clock, bcs, fields)z
end

function initialize!(state, model::VegetationSoilModel)
    initialize!(state, model, model.initializer)
    initialize!(state, model.surface_energy_balance)
    initialize!(state, model.surface_hydrology)
    initialize!(state, model.vegetation)
    initialize!(state, model.soil)
end

function compute_auxiliary!(state, model::VegetationSoilModel)
    atmos = get_atmosphere(model)
    seb = get_surface_energy_balance(model)
    surface_hydrology = get_surface_hydrology(model)
    compute_auxiliary!(state, model, atmos)
    compute_auxiliary!(state, model.soil)
    compute_auxiliary(state, seb)
    compute_auxiliary!(state, model, model.vegetation)
    compute_auxiliary!(state, model, surface_hydrology)
    # recompute surface energy fluxes
    compute_surface_energy_fluxes!(state, model, seb)
end

function compute_tendencies!(state, model::VegetationSoilModel)
    compute_tendencies!(state, model, model.atmosphere)
    compute_tendencies!(state, model.surface_energy_balance)
    compute_tendencies!(state, model, model.surface_hydrology)
    compute_tendencies!(state, model.soil)
    compute_tendencies!(state, model, model.vegetation)
end
