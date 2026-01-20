@kwdef struct VegetationSoilModel{
    NF,
    GridType<:AbstractLandGrid{NF},
    Atmosphere<:AbstractAtmosphere,
    SEB<:AbstractSurfaceEnergyBalance,
    Vegetation<:AbstractVegetationModel{NF, GridType},
    SoilModel<:AbstractSoilModel{NF, GridType},
    Initializer<:AbstractInitializer,
} <: AbstractLandModel{NF, GridType}
    "Spatial discretization"
    grid::GridType

    "Near-surface atmospheric conditions"
    atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Surface energy balance"
    surface_energy::SEB = SurfaceEnergyBalance(eltype(grid))

    "Vegetation model"
    vegetation::Vegetation = VegetationModel(grid)

    "Soil model"
    soil::SoilModel = SoilModel(grid)

    "Initializer for coupled model"
    initializer::Initializer = DefaultInitializer()
end

get_surface_energy_balance(model::VegetationSoilModel) = model.surface_energy

get_soil_energy_balance(model::VegetationSoilModel) = model.soil.energy

get_soil_hydrology(model::VegetationSoilModel) = model.soil.hydrology

get_soil_stratigraphy(model::VegetationSoilModel) = model.soil.strat

get_soil_biogeochemistry(model::VegetationSoilModel) = model.soil.biogeochem

# just ignore that these are duplicated in veg model for now...
get_constants(model::VegetationSoilModel) = model.soil.constants

function variables(model::VegetationSoilModel)
    atmos_vars = variables(model.atmosphere)
    veg_vars = variables(model.vegetation)
    soil_vars = variables(model.soil)
    return tuplejoin(atmos_vars, veg_vars, soil_vars)
end

function initialize(
    model::VegetationSoilModel{NF};
    clock = Clock(time=zero(NF)),
    boundary_conditions = (;),
    fields = (;),
    external_variables = ()
) where {NF}
    vars = Variables(variables(model)..., external_variables...)
    # TODO
    return initialize(vars, model.grid, clock)
end

function initialize!(state, model::VegetationSoilModel)
    initialize!(state, model, model.initializer)
    initialize!(state, model.soil)
    initialize!(state, model.soil)
end

function compute_auxiliary!(state, model::VegetationSoilModel)
    grid = get_grid(model)
    compute_auxiliary!(state, model, model.atmosphere)
    compute_auxiliary!(state, model, model.vegetation)
    compute_auxiliary!(state, model.soil)
end

function compute_tendencies!(state, model::VegetationSoilModel)
    compute_tendencies!(state, model, model.vegetation)
    compute_tendencies!(state, model.soil)
end
