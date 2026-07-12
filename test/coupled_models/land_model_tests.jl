using Terrarium
using Test

using Oceananigans.BoundaryConditions: BoundaryCondition, Flux

@testset "LandModel: Soil, no vegetation" begin
    grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 50))
    swrc = VanGenuchten(α = 2.0, n = 2.0)
    hydraulic_properties = ConstantSoilHydraulics(eltype(grid); swrc, unsat_hydraulic_cond = UnsatKVanGenuchten(eltype(grid)))
    hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties)
    soil = SoilEnergyWaterCarbon(eltype(grid); hydrology)
    land = LandModel(grid; soil, vegetation = nothing)
    # Test vegetation = nothing results in bare ground evaporation and no canopy interception
    @test land.surface_hydrology.evapotranspiration isa Terrarium.BareGroundEvaporation
    @test land.surface_hydrology.canopy_interception isa Terrarium.NoCanopyInterception
    # Variably saturated with water table
    initializers = (
        temperature = (x, z) -> 5.0 - 0.02 * z,
        saturation_water_ice = (x, z) -> min(1, 0.8 - 0.05 * z),
    )
    integrator = initialize(land; initializers)
    # Check that infiltration is correctly coupled to soil hydrology
    set!(integrator.state.infiltration, 1.0e-8)
    sat_top_bc = integrator.state.saturation_water_ice.boundary_conditions.top
    @test isa(sat_top_bc, BoundaryCondition{<:Flux})
    @test all(Field(sat_top_bc.condition) .≈ -1.0e-8) # note the negative sign
    # Check that ground heat flux is coupled to soil energy
    energy_top_bc = integrator.state.internal_energy.boundary_conditions.top
    @test isa(energy_top_bc, BoundaryCondition{<:Flux})
    @test energy_top_bc.condition == integrator.state.ground_heat_flux
    # Advance one timestep
    timestep!(integrator, 60.0)
    @test all(isfinite.(integrator.state.saturation_water_ice))
    @test all(isfinite.(integrator.state.internal_energy))
    @test all(isfinite.(integrator.state.ground_heat_flux))
end

@testset "LandModel: Coupled vegetation-soil" begin
    grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 50))
    swrc = VanGenuchten(α = 2.0, n = 2.0)
    hydraulic_properties = ConstantSoilHydraulics(eltype(grid); swrc, unsat_hydraulic_cond = UnsatKVanGenuchten(eltype(grid)))
    hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties)
    soil = SoilEnergyWaterCarbon(eltype(grid); hydrology)
    vegetation = VegetationCarbon(eltype(grid))
    land = LandModel(grid; soil, vegetation)
    @test land.surface_hydrology.evapotranspiration isa Terrarium.PALADYNCanopyEvapotranspiration
    @test land.surface_hydrology.canopy_interception isa Terrarium.PALADYNCanopyInterception
    # Variably saturated with water table
    initializers = (
        temperature = (x, z) -> 5.0 - 0.02 * z,
        saturation_water_ice = (x, z) -> min(1, 0.8 - 0.05 * z),
        carbon_vegetation = 0.1,
    )
    integrator = initialize(land; initializers)
    # Check that infiltration is correctly coupled to soil hydrology
    set!(integrator.state.infiltration, 1.0e-8)
    sat_top_bc = integrator.state.saturation_water_ice.boundary_conditions.top
    @test isa(sat_top_bc, BoundaryCondition{<:Flux})
    @test all(Field(sat_top_bc.condition) .≈ -1.0e-8) # note the negative sign
    # Check that ground heat flux is coupled to soil energy
    energy_top_bc = integrator.state.internal_energy.boundary_conditions.top
    @test isa(energy_top_bc, BoundaryCondition{<:Flux})
    @test energy_top_bc.condition == integrator.state.ground_heat_flux
    # Advance one timestep
    timestep!(integrator, 60.0)
    @test all(isfinite.(integrator.state.saturation_water_ice))
    @test all(isfinite.(integrator.state.internal_energy))
    @test all(isfinite.(integrator.state.ground_heat_flux))
    @test all(isfinite.(integrator.state.carbon_vegetation))
    # TODO: also check ET and veg processes once they are working...
end

@testset "LandModel: Snow-coupled soil" begin
    grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 50))
    swrc = VanGenuchten(α = 2.0, n = 2.0)
    hydraulic_properties = ConstantSoilHydraulics(eltype(grid); swrc, unsat_hydraulic_cond = UnsatKVanGenuchten(eltype(grid)))
    hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties)
    soil = SoilEnergyWaterCarbon(eltype(grid); hydrology)
    land = LandModel(grid; soil, snow = SingleLayerSnow(eltype(grid)), vegetation = nothing)
    @test land.snow isa SingleLayerSnow
    # Frozen soil beneath a cold snowpack
    initializers = (
        temperature = (x, z) -> -1.0 - 0.02 * z,
        saturation_water_ice = (x, z) -> min(1, 0.8 - 0.05 * z),
        snow_water_equivalent = 0.2,
        snow_temperature = -5.0,
    )
    integrator = initialize(land; initializers)
    # With snow, the soil-top energy BC reads the blended `soil_heat_flux`, not `ground_heat_flux` directly
    energy_top_bc = integrator.state.internal_energy.boundary_conditions.top
    @test isa(energy_top_bc, BoundaryCondition{<:Flux})
    @test energy_top_bc.condition === integrator.state.soil_heat_flux
    # Advance one timestep
    timestep!(integrator, 60.0)
    # A 0.2 m SWE pack should be essentially fully snow-covered (cover fraction is diagnosed in compute_auxiliary!)
    @test all(integrator.state.snow_cover_fraction .> 0.9)
    @test all(isfinite.(integrator.state.snow_energy))
    @test all(isfinite.(integrator.state.snow_water_equivalent))
    @test all(isfinite.(integrator.state.snow_temperature))
    @test all(isfinite.(integrator.state.basal_heat_flux))
    @test all(isfinite.(integrator.state.soil_heat_flux))
    @test all(isfinite.(integrator.state.internal_energy))
end
