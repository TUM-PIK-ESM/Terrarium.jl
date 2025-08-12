using Terrarium
using Test

using Terrarium: getproperties
using SpecialFunctions: erfc

import Oceananigans.BoundaryConditions: ValueBoundaryCondition, FluxBoundaryCondition, NoFluxBoundaryCondition, regularize_field_boundary_conditions

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
    # model setup
    grid = ColumnGrid(ExponentialSpacing(Δz_min=0.05, Δz_max=1.0, N=50))
    # specify periodic upper boundary condition
    upperbc(x, y, t) = sin(2π*t / P)
    boundary_conditions = VarBoundaryConditions(:temperature, top=ValueBoundaryCondition(upperbc))
    # temperature initial condition
    initializer = VarInitializer(T₀, :temperature)
    model = SoilModel(; grid, initializer, boundary_conditions)
    sim = initialize(model)
    # timestep!(sim)
    # TODO add check
end

@testset "Step heat diffusion" begin
    # Temperatures should be > 0 to avoid effects from phase change
    T₀ = 1.0
    T₁ = 2.0
    # model setup
    grid = ColumnGrid(UniformSpacing(Δz=0.01, N=50))
    z_centers = grid.grid.z.cᵃᵃᶜ[1:end-3]
    # temperature initial condition
    initializer = VarInitializer(T₀, :temperature)
    # constant upper boundary temperature set to T₁
    boundary_conditions = VarBoundaryConditions(:temperature, top=ValueBoundaryCondition(T₁))
    # set carbon content to zero so the soil has only a mineral constituent
    biogeochem = ConstantSoilCarbonDenisty(ρ_soc=0.0)
    # set porosity to zero to remove influence of pore space;
    # this is just a hack to configure the model to simulate heat conduction in a fully solid medium
    hydraulic_properties = PrescribedHydraulics(porosity=0.0)
    hydrology = SoilHydrology(; hydraulic_properties)
    model = SoilModel(; grid, hydrology, biogeochem, initializer, boundary_conditions)
    sim = initialize(model)
    # TODO: Rewrite this part once we have a proper output handling system
    Ts_buf = [deepcopy(sim.state.temperature)]
    ts = [0.0]
    dt = 1.0
    # run for one hour, saving every time step
    while current_time(sim) < 3600
        timestep!(sim, dt)
        push!(Ts_buf, deepcopy(sim.state.temperature))
        push!(ts, current_time(sim))
    end

    soil_thermal_props = model.energy.thermal_properties
    k = soil_thermal_props.cond.mineral
    c = soil_thermal_props.heatcap.mineral
    α = k / c
    ΔT_sol = heat_conduction_linear_step_ub(T₁ - T₀, α)
    Ts = reduce(hcat, Ts_buf)[1,:,:]
    Ts_target = T₀ .+ ΔT_sol.(reshape(abs.(z_centers), 1, :), reshape(ts, :, 1))
    error = Ts .- Ts_target

    # Check error at last time step
    last_error_threshold = 1e-3
    @test maximum(abs.(error[end,:])) < last_error_threshold

    # Check total error over all time steps
    # TODO: This threshold is too high, need to double check this and see if there is an actual bug.
    total_error_threshold = 0.5
    @test maximum(abs.(error)) < total_error_threshold
end

@testset "Thermal properties" begin
    thermal_props = SoilThermalProperties()
    
    # check that all necessary properties are defined
    @test all(map(∈(propertynames(thermal_props.cond)), (:water, :ice, :air, :mineral, :organic)))
    @test all(map(∈(propertynames(thermal_props.heatcap)), (:water, :ice, :air, :mineral, :organic)))
    
    # check that all default values are valid
    @test all(>(0), getproperties(thermal_props.cond))
    @test all(>(0), getproperties(thermal_props.heatcap))
    
    # sanity checks for bulk thermal properties
    @test Terrarium.thermalconductivity(thermal_props, (water=1.0, ice=0.0, air=0.0, mineral=0.0, organic=0.0)) ≈ thermal_props.cond.water
    @test Terrarium.thermalconductivity(thermal_props, (water=0.0, ice=1.0, air=0.0, mineral=0.0, organic=0.0)) ≈ thermal_props.cond.ice
    @test Terrarium.thermalconductivity(thermal_props, (water=0.0, ice=0.0, air=1.0, mineral=0.0, organic=0.0)) ≈ thermal_props.cond.air
    @test Terrarium.thermalconductivity(thermal_props, (water=0.0, ice=0.0, air=0.0, mineral=1.0, organic=0.0)) ≈ thermal_props.cond.mineral
    @test Terrarium.thermalconductivity(thermal_props, (water=0.0, ice=0.0, air=0.0, mineral=0.0, organic=1.0)) ≈ thermal_props.cond.organic
end
