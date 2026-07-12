"""
    $TYPEDEF

Fully-coupled land model integrating atmosphere, surface energy balance, surface hydrology,
vegetation, and soil processes.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct LandModel{
        NF,
        GridType <: AbstractLandGrid{NF},
        Vegetation <: Optional{AbstractVegetation{NF}},
        Soil <: AbstractSoil{NF},
        Snow <: Optional{AbstractSnow{NF}},
        SEB <: AbstractSurfaceEnergyBalance,
        Hydrology <: AbstractSurfaceHydrology,
        Atmosphere <: AbstractAtmosphere,
        Initializer <: AbstractInitializer,
        Timestepper <: AbstractTimeStepper{NF},
    } <: AbstractLandModel{NF, GridType}
    "Spatial discretization"
    grid::GridType

    "Vegetation processes"
    @component vegetation::Vegetation = VegetationCarbon(eltype(grid))

    "Soil processes"
    @component soil::Soil = default_soil(grid, vegetation)

    "Snow processes (default `nothing`: no snowpack)"
    @component snow::Snow = nothing

    "Surface energy balance"
    @component surface_energy_balance::SEB = default_surface_energy_balance(grid, vegetation, soil)

    "Surface hydrology scheme"
    @component surface_hydrology::Hydrology = default_surface_hydrology(grid, vegetation, soil)

    "Near-surface atmospheric conditions"
    @component atmosphere::Atmosphere = PrescribedAtmosphere(eltype(grid))

    "Physical constants"
    @component constants::PhysicalConstants{NF} = PhysicalConstants(eltype(grid))

    "State variable initializer"
    @component initializer::Initializer = DefaultInitializer(eltype(grid))

    "Time stepper: a single `AbstractTimeStepper` (e.g. `ForwardEuler`, `Heun`) or an `IMEX`"
    @component timestepper::Timestepper = default_timestepper(eltype(grid))
end

function StateVariables(
        model::LandModel{NF};
        clock = Clock(time = zero(NF)),
        boundary_conditions = (;),
        fields = (;),
        input_variables = ()
    ) where {NF}
    grid = get_grid(model)
    interface_vars = interface_variables(model)
    # Snow adds a blended `soil_heat_flux` coupling auxiliary to the state when present
    vars = Variables(variables(model)..., interface_vars..., input_variables...)
    # Initialize BC fields for coupling
    ground_heat_flux = initialize(vars.auxiliary.ground_heat_flux, grid, clock, boundary_conditions, fields)
    infiltration = initialize(vars.auxiliary.infiltration, grid, clock, boundary_conditions, fields)
    # Initialize soil heat flux
    soil_heat_flux = if isnothing(model.snow)
        ground_heat_flux
    else
        initialize(vars.auxiliary.soil_heat_flux, grid, clock, boundary_conditions, fields)
    end
    soil_heat_flux_bc = SoilHeatFlux(soil_heat_flux)
    # Note that the hydrology module computes infiltration as positive so we need to negate it here
    # since fluxes are by convention positive upwards
    infiltration_bc = InfiltrationFlux(-infiltration)
    bcs = merge_boundary_conditions(boundary_conditions, soil_heat_flux_bc, infiltration_bc)
    # The snow energy tendency reads its boundary heat fluxes from the `surface_heat_flux`/`basal_heat_flux`
    # input fields; alias them to the surface energy balance closure flux (`ground_heat_flux`) and the
    # blended soil-top flux (`soil_heat_flux`) so no additional state fields are allocated (no-op without snow).
    snow_flux_aliases = isnothing(model.snow) ? (;) : (; surface_heat_flux = ground_heat_flux, basal_heat_flux = soil_heat_flux)
    # Merge user-defined fields with BC fields
    fields = merge((; ground_heat_flux, infiltration, soil_heat_flux), snow_flux_aliases, fields)
    return StateVariables(vars, grid; clock, timestepper = get_timestepper(model), model, boundary_conditions = bcs, fields)
end

interface_variables(::LandModel) = (
    auxiliary(:soil_heat_flux, XY(); units = u"W/m^2", desc = "Blended heat flux into the soil top (snow base + bare ground)"),
)

function initialize!(state, model::LandModel)
    initialize!(state, model, model.initializer)
    grid = get_grid(model)
    initialize!(state, grid, model.surface_hydrology)
    # TODO: change when refactoring model/process types
    initialize!(state, grid, model.vegetation, model.constants, model.atmosphere)
    initialize!(state, grid, model.soil, model.constants)
    # Initialize the snow (invclosure: T, W -> E) after the soil; no-op without snow (`nothing`)
    initialize!(state, grid, model.snow, model.constants)
    # Initialize the SEB after the soil so that ground_temperature is available
    initialize!(state, grid, model.surface_energy_balance)
    return nothing
end

function compute_auxiliary!(state, model::LandModel)
    grid = get_grid(model)
    compute_auxiliary!(state, grid, model.atmosphere)
    compute_auxiliary!(state, grid, model.soil, model.constants)
    compute_auxiliary!(state, grid, model.snow, model.constants)
    compute_auxiliary!(state, grid, model.vegetation, model.constants, model.atmosphere, model.soil)
    compute_auxiliary!(state, grid, model.surface_hydrology, model.constants, model.atmosphere, model.soil, model.vegetation)
    compute_auxiliary!(state, grid, model.surface_energy_balance, model.constants, model.atmosphere, model.surface_hydrology, model.snow)
    # Blend the soil-top heat flux from the snow base and bare ground after the SEB; no-op without snow
    compute_soil_snow_fluxes!(state, grid, model.snow, model.constants)
    return nothing
end

function compute_tendencies!(state, model::LandModel)
    grid = get_grid(model)
    compute_tendencies!(state, grid, model.surface_hydrology)
    compute_tendencies!(state, grid, model.soil, model.constants)
    compute_tendencies!(state, grid, model.vegetation)
    # Snow tendencies after surface hydrology; no-op without snow (`nothing`)
    compute_tendencies!(state, grid, model.snow, model.constants, model.atmosphere)
    return nothing
end

function closure!(state, model::LandModel)
    grid = get_grid(model)
    closure!(state, grid, model.soil, model.constants)
    closure!(state, grid, model.snow, model.constants)
    return nothing
end

function invclosure!(state, model::LandModel)
    grid = get_grid(model)
    invclosure!(state, grid, model.soil, model.constants)
    invclosure!(state, grid, model.snow, model.constants)
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
