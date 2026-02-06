using Terrarium
using Test

using Terrarium: getproperties
using SpecialFunctions: erfc

import Oceananigans.BoundaryConditions: ValueBoundaryCondition, FluxBoundaryCondition, NoFluxBoundaryCondition, regularize_field_boundary_conditions

@testset "Thermal properties" begin
    thermal_props = SoilThermalProperties(Float64)
    
    # check that all necessary properties are defined
    @test all(map(∈(propertynames(thermal_props.conductivities)), (:water, :ice, :air, :mineral, :organic)))
    @test all(map(∈(propertynames(thermal_props.heat_capacities)), (:water, :ice, :air, :mineral, :organic)))
    
    # check that all default values are valid
    @test all(>(0), getproperties(thermal_props.conductivities))
    @test all(>(0), getproperties(thermal_props.heat_capacities))
    
    # sanity checks for bulk thermal properties
    @test Terrarium.compute_thermal_conductivity(thermal_props, SoilVolume(porosity=1.0, saturation=1.0, liquid=1.0)) ≈ thermal_props.conductivities.water
    @test Terrarium.compute_thermal_conductivity(thermal_props, SoilVolume(porosity=1.0, saturation=1.0, liquid=0.0)) ≈ thermal_props.conductivities.ice
    @test Terrarium.compute_thermal_conductivity(thermal_props, SoilVolume(porosity=1.0, saturation=0.0, liquid=0.0)) ≈ thermal_props.conductivities.air
    @test Terrarium.compute_thermal_conductivity(thermal_props, SoilVolume(porosity=0.0, saturation=0.0)) ≈ thermal_props.conductivities.mineral
    @test Terrarium.compute_thermal_conductivity(thermal_props, SoilVolume(porosity=0.0, saturation=0.0, solid = MineralOrganic(organic = 1.0))) ≈ thermal_props.conductivities.organic
end

@testset "Soil energy: initialize!" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing())
    energy = SoilEnergyBalance(eltype(grid))
    soil = SoilEnergyWaterCarbon(eltype(grid); energy)
    constants = PhysicalConstants(eltype(grid))
    state = initialize(soil, grid)
    # Test initialization with 0°C temperature
    set!(state.temperature, 0.0)
    Terrarium.initialize!(state, grid, energy, soil, constants)
    @test all(state.liquid_water_fraction .≈ 1.0)
    @test all(state.internal_energy .≈ 0.0)
    # Test initialization with positive temperature
    set!(state.temperature, 1.0) # initialize temperature to zero
    Terrarium.initialize!(state, grid, energy, soil, constants)
    @test all(state.liquid_water_fraction .≈ 1.0)
    @test all(state.internal_energy .> 0.0)
    # Test initialization with negative temperature
    set!(state.temperature, -1.0) # initialize temperature to zero
    Terrarium.initialize!(state, grid, energy, soil, constants)
    @test all(state.liquid_water_fraction .≈ 0.0)
    @test all(state.internal_energy .< 0.0)
end

@testset "Soil energy: compute_tendencies!" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))
    energy = SoilEnergyBalance(eltype(grid))
    soil = SoilEnergyWaterCarbon(eltype(grid); energy)
    constants = PhysicalConstants(eltype(grid))
    state = initialize(soil, grid)
    set!(state.temperature, (x, z) -> 0.0 - 0.01 * z)
    Terrarium.initialize!(state, grid, energy, soil, constants)
    compute_tendencies!(state, grid, energy, soil)
    @test all(isfinite.(state.tendencies.internal_energy))
end

@testset "Soil energy: closure!" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))
    energy = SoilEnergyBalance(eltype(grid))
    soil = SoilEnergyWaterCarbon(eltype(grid); energy)
    constants = PhysicalConstants(eltype(grid))
    state = initialize(soil, grid)
    set!(state.internal_energy, 1e6)
    Terrarium.closure!(state, grid, energy, soil, constants)
    @test all(state.temperature .> 0)
    @test all(state.liquid_water_fraction .≈ 1.0)
end

# Analytical solution to the 1D heat equation with diffusivity α
# and periodic upper boundary with mean T₀, amplitude A, and period P.
function heat_conduction_linear_periodic_ub(T₀, A, P, α)
    T(z,t) = T₀ + A*exp(-z*sqrt(π/(α*P)))*sin(2π*t/P - z*sqrt(π/(α*P)))
    return T
end

# Analytical solution to the 1D heat equation with diffusivity α
# under a step change of the upper boundary condition at steady state.
function heat_conduction_linear_step_ub(ΔTₛ, α)
    ΔT(z,t) = ΔTₛ * erfc(z / (2*sqrt(α*t)))
    return ΔT
end

@testset "Heat diffusion with periodic upper BC" begin
    # parameters
    T₀ = 2.0
    A = 1.0
    P = 24*3600
    k = 2.0
    c = 1e6
    α = k / c
    T_sol = heat_conduction_linear_periodic_ub(T₀, A, P, α)

    # model setup
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(Δz_min=0.05, Δz_max=100.0, N=100))
    # temperature initial condition
    initializer = FieldInitializers(
        temperature = (x,z) -> T_sol(-z, 0.0)
    )
    # set carbon content to zero so the soil has only a mineral constituent
    biogeochem = ConstantSoilCarbonDensity(ρ_soc=0.0)
    # set porosity to zero to remove influence of pore space;
    # this is just a hack to configure the model to simulate heat conduction in a fully solid medium
    soil_porosity = ConstantSoilPorosity(mineral_porosity=0.0)
    strat = HomogeneousStratigraphy(Float64; porosity = soil_porosity)
    # set thermal properties
    thermal_properties = SoilThermalProperties(
        eltype(grid);
        conductivities=SoilThermalConductivities(mineral=k),
        heat_capacities=SoilHeatCapacities(mineral=c),
    )
    energy = SoilEnergyBalance(eltype(grid); thermal_properties)
    soil = SoilEnergyWaterCarbon(eltype(grid); energy, strat, biogeochem)
    model = SoilModel(grid; soil, initializer)
    # periodic upper boundary temperature
    upperbc(z, t) = T₀ + A*sin(2π*t/P)
    bcs = PrescribedSurfaceTemperature(:Tsurf, upperbc)
    integrator = initialize(model, ForwardEuler(), boundary_conditions = bcs)
    # TODO: Rewrite this part once we have a proper output handling system
    Ts_buf = [deepcopy(integrator.state.temperature)]
    ts = [0.0]
    Δt = 60.0
    # run for one hour, saving every time step
    while current_time(integrator) < 2*P
        timestep!(integrator, Δt)
        push!(Ts_buf, deepcopy(integrator.state.temperature))
        push!(ts, current_time(integrator))
    end

    z_centers = znodes(integrator.state.temperature)
    Ts = reduce(hcat, Ts_buf)[1,:,:]
    Ts_target = T_sol.(reshape(-z_centers, 1, :), reshape(ts, :, 1))
    relative_error = abs.((Ts .- Ts_target) ./ Ts_target)
    @show maximum(relative_error)
    @test maximum(relative_error) < 0.1
end

@testset "Step heat diffusion" begin
    # Temperatures should be > 0 to avoid effects from phase change
    T₀ = 1.0
    T₁ = 2.0
    # model setup
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(Δz_min=0.01, Δz_max=100.0, N=100))
    # temperature initial condition
    initializer = FieldInitializers(temperature=T₀)
    # set carbon content to zero so the soil has only a mineral constituent
    biogeochem = ConstantSoilCarbonDensity(ρ_soc=0.0)
    # set porosity to zero to remove influence of pore space;
    # this is just a hack to configure the model to simulate heat conduction in a fully solid medium
    soil_porosity = ConstantSoilPorosity(mineral_porosity=0.0)
    strat = HomogeneousStratigraphy(Float64; porosity = soil_porosity)
    soil = SoilEnergyWaterCarbon(eltype(grid); strat, biogeochem)
    model = SoilModel(grid; soil, initializer)
    # constant upper boundary temperature set to T₁
    bcs = (temperature = (top = ValueBoundaryCondition(T₁),),)
    integrator = initialize(model, ForwardEuler(), boundary_conditions = bcs)
    # TODO: Rewrite this part once we have a proper output handling system
    Ts_buf = [deepcopy(integrator.state.temperature)]
    ts = [0.0]
    Δt = 10.0
    # run for 24 hours, saving every time step
    while current_time(integrator) < 24*3600
        timestep!(integrator, Δt)
        push!(Ts_buf, deepcopy(integrator.state.temperature))
        push!(ts, current_time(integrator))
    end

    soil_thermal_props = soil.energy.thermal_properties
    k = soil_thermal_props.conductivities.mineral
    c = soil_thermal_props.heat_capacities.mineral
    α = k / c
    z_centers = znodes(integrator.state.temperature)
    ΔT_sol = heat_conduction_linear_step_ub(T₁ - T₀, α)
    Ts = reduce(hcat, Ts_buf)[1,:,:]
    Ts_target = T₀ .+ ΔT_sol.(reshape(-z_centers, 1, :), reshape(ts, :, 1))
    relative_error = abs.((Ts .- Ts_target) ./ Ts_target)
    @show maximum(relative_error)

    # Check error at last time step
    last_error_threshold = 1e-3
    @test maximum(relative_error[end,:]) < last_error_threshold

    # Check total error over all time steps
    max_error_threshold = 0.1
    @test maximum(relative_error) < max_error_threshold
end
