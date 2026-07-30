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
    # Cover fraction is diagnosed in compute_auxiliary! as f = W/(W + W_ref); with the default
    # half-coverage W_ref = 0.1 m, a 0.2 m SWE pack is majority-covered (f = 0.2/0.3 ≈ 0.67)
    @test all(integrator.state.snow_cover_fraction .> 0.6)
    @test all(isfinite.(integrator.state.snow_energy))
    @test all(isfinite.(integrator.state.snow_water_equivalent))
    @test all(isfinite.(integrator.state.snow_temperature))
    @test all(isfinite.(integrator.state.basal_heat_flux))
    @test all(isfinite.(integrator.state.soil_heat_flux))
    @test all(isfinite.(integrator.state.internal_energy))
end

@testset "LandModel: snow meltwater reaches the soil surface" begin
    # Water-conservation check for the snow→soil coupling (Option A): with a melting snowpack and no
    # rain, the meltwater outflow M_r is the sole surface water input, so the runoff scheme must
    # partition exactly M_r into infiltration + surface runoff (no drainage when the surface reservoir
    # is empty): infiltration + surface_runoff ≈ M_r.
    grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 50))
    soil = SoilEnergyWaterCarbon(eltype(grid); hydrology = SoilHydrology(eltype(grid), RichardsEq()))
    snow = SingleLayerSnow(eltype(grid))
    land = LandModel(grid; soil, snow, vegetation = nothing)
    # Fully-liquid (melting) pack over unsaturated soil; no rainfall
    initializers = (
        temperature = (x, z) -> 1.0 - 0.02 * z,
        saturation_water_ice = (x, z) -> 0.5,
        snow_water_equivalent = 0.2,
        snow_temperature = 0.0,   # invclosure -> fully-liquid pack (θ_liq = 1)
    )
    integrator = initialize(land; initializers)
    state = integrator.state
    set!(state.rainfall, 0.0)
    Terrarium.closure!(state, land)
    compute_auxiliary!(state, land)
    # meltwater outflow implied by the diagnosed liquid fraction
    θ_liq = Array(interior(state.snow_liquid_fraction))
    M_r = Terrarium.compute_meltwater_outflow.(Ref(snow.hydraulic_properties), θ_liq)
    @test all(θ_liq .≈ 1)                    # melting pack
    @test all(M_r .> 0)                       # meltwater draining
    infil = Array(interior(state.infiltration))
    runoff = Array(interior(state.surface_runoff))
    @test all(isapprox.(infil .+ runoff, M_r; rtol = 1.0e-6))
end

@testset "LandModel: latent flux partitioned between ground evaporation and sublimation" begin
    # The surface latent flux is split by snow-covered fraction: ground evaporation over (1 − f_snow),
    # snow sublimation over f_snow. The ground evaporation is scaled at its source by (1 − f_snow), so the
    # stored `evaporation_ground` must equal (1 − f_snow)·g·Δq(T_skin) — no double-counting with sublimation.
    grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 50))
    soil = SoilEnergyWaterCarbon(eltype(grid); hydrology = SoilHydrology(eltype(grid), RichardsEq()))
    inits = (
        temperature = (x, z) -> 2.0 - 0.02 * z, saturation_water_ice = (x, z) -> 0.8,
        snow_water_equivalent = 0.5, snow_temperature = -2.0,
    )
    land = LandModel(grid; soil, snow = SingleLayerSnow(eltype(grid)), vegetation = nothing)
    it = initialize(land; initializers = inits)
    Terrarium.closure!(it.state, land)
    compute_auxiliary!(it.state, land)
    thermo = land.constants.thermodynamics
    f = Array(interior(it.state.snow_cover_fraction))
    E = Array(interior(it.state.evaporation_ground))
    g = Array(interior(it.state.ground_evaporation_conductance))
    Ts = Array(interior(it.state.skin_temperature))
    p = Array(interior(it.state.air_pressure))
    q_air = Array(interior(it.state.specific_humidity))
    # unscaled bulk-aerodynamic ground evaporation g·Δq(T_skin); the stored flux is this scaled by (1 − f)
    Δq = Terrarium.specific_humidity_difference.(Ref(thermo), p, q_air, Ts)
    @test all(f .> 0.8)                                      # 0.5 m SWE -> f = 0.5/(0.5 + 0.1) ≈ 0.83 (W_ref = 0.1 m)
    @test all(isapprox.(E, (1 .- f) .* g .* Δq; rtol = 1.0e-6))   # ground evaporation scaled by (1 − f_snow)
    @test all(isfinite.(interior(it.state.sublimation)))
    @test all(isfinite.(interior(it.state.latent_heat_flux)))
end
