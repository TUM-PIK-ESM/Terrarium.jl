"""
    $TYPEDEF

Fully-coupled land model integrating atmosphere, surface energy balance, surface hydrology,
vegetation, and soil processes.

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct LandModel{
        NF,
        GridType <: AbstractLandGrid{NF},
        Vegetation <: Optional{AbstractVegetation{NF}},
        Soil <: AbstractSoil{NF},
        SEB <: AbstractSurfaceEnergyBalance,
        Hydrology <: AbstractSurfaceHydrology,
        Atmosphere <: AbstractAtmosphere,
        Constants <: PhysicalConstants{NF},
        Initializer <: AbstractInitializer,
    } <: AbstractLandModel{NF, GridType}
    "Spatial discretization"
    grid::GridType

    "Vegetation processes"
    vegetation::Vegetation = VegetationCarbon(eltype(grid))

    "Soil processes"
    soil::Soil = default_soil(grid, vegetation)

    "Surface energy balance"
    surface_energy_balance::SEB = default_surface_energy_balance(grid, vegetation, soil)

    "Surface hydrology scheme"
    surface_hydrology::Hydrology = default_surface_hydrology(grid, vegetation, soil)

    "Near-surface atmospheric conditions"
    atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Physical constants"
    constants::Constants = PhysicalConstants(eltype(grid))

    "State variable initializer"
    initializer::Initializer = DefaultInitializer(eltype(grid))
end

function initialize(
        model::LandModel{NF};
        clock = Clock(time = zero(NF)),
        boundary_conditions = (;),
        fields = (;),
        input_variables = ()
    ) where {NF}
    grid = get_grid(model)
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
    return initialize(vars, grid; clock, boundary_conditions = bcs, fields)
end

function initialize!(state, model::LandModel)
    initialize!(state, model, model.initializer)
    grid = get_grid(model)
    initialize!(state, grid, model.surface_energy_balance)
    initialize!(state, grid, model.surface_hydrology)
    # TODO: change when refactoring model/process types
    initialize!(state, grid, model.vegetation, model.atmosphere, model.constants)
    initialize!(state, grid, model.soil, model.constants)
    return nothing
end

function compute_auxiliary!(state, model::LandModel)
    grid = get_grid(model)
    compute_auxiliary!(state, grid, model.atmosphere)
    compute_auxiliary!(state, grid, model.soil, model.constants)
    compute_auxiliary!(state, grid, model.vegetation, model.atmosphere, model.constants, model.soil)
    compute_auxiliary!(state, grid, model.surface_hydrology, model.atmosphere, model.constants, model.vegetation, model.soil)
    compute_auxiliary!(state, grid, model.surface_energy_balance, model.atmosphere, model.constants, model.surface_hydrology)
    compute_surface_energy_fluxes!(state, grid, model.surface_energy_balance, model.atmosphere, model.constants, model.surface_hydrology)
    return nothing
end

function compute_tendencies!(state, model::LandModel)
    grid = get_grid(model)
    compute_tendencies!(state, grid, model.surface_hydrology)
    compute_tendencies!(state, grid, model.soil, model.constants)
    compute_tendencies!(state, grid, model.vegetation)
    return nothing
end

function closure!(state, model::LandModel)
    grid = get_grid(model)
    closure!(state, grid, model.soil, model.constants)
    return nothing
end

function invclosure!(state, model::LandModel)
    grid = get_grid(model)
    invclosure!(state, grid, model.soil, model.constants)
    return nothing
end

# Default soil types
default_soil(grid::AbstractLandGrid, ::Nothing) = SoilEnergyWaterCarbon(eltype(grid))
default_soil(grid::AbstractLandGrid, ::AbstractVegetation) = SoilEnergyWaterCarbon(eltype(grid), hydrology = SoilHydrology(eltype(grid), RichardsEq()))

# Default SEB
default_surface_energy_balance(grid::AbstractLandGrid, vegetation, soil) = SurfaceEnergyBalance(eltype(grid))

# Default surface hydrology
default_surface_hydrology(grid::AbstractLandGrid, ::AbstractVegetation, soil) = SurfaceHydrology(eltype(grid))
function default_surface_hydrology(grid::AbstractLandGrid, ::Nothing, soil)
    return SurfaceHydrology(
        eltype(grid),
        evapotranspiration = BareGroundEvaporation(eltype(grid)),
        canopy_interception = NoCanopyInterception(eltype(grid))
    )
end
