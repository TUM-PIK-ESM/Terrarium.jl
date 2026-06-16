using Terrarium
using Thermodynamics
using Oceananigans
using Test

@testset "Prescribed skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))
    skin_temperature = PrescribedSkinTemperature(eltype(grid))
    seb = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, surface_energy_balance = seb)
    state = initialize(model)
    @test hasproperty(state.inputs, :skin_temperature)
    set!(state.skin_temperature, 1.0)
    compute_auxiliary!(state, grid, skin_temperature)
    @test all(state.skin_temperature .≈ 1.0)
end

@testset "Implicit skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))
    clock = Clock(time = 0.0)
    skin_temperature = ImplicitSkinTemperature(κₛ = 0.5)
    seb = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, surface_energy_balance = seb)
    state = initialize(model)
    @test !hasproperty(state.inputs, :skin_temperature)
    @test hasproperty(state.inputs, :ground_temperature)

    # Sunny and dry in Bergen Norway (Figure 5.11, Shuttleworth 2012)
    set!(state.surface_shortwave_down, 600.0) # sunny conditions
    set!(state.surface_longwave_down, 300.0)
    air_pressure = 101_325 # standard pressure
    set!(state.air_pressure, air_pressure) # standard pressure
    air_temperature = 10.0 # 10 °C
    set!(state.air_temperature, air_temperature)
    thermodyn_constants = model.constants.thermodynamics
    air_density = Thermodynamics.air_density(thermodyn_constants, air_temperature + 273.15, air_pressure)
    specific_humidity_saturation = Terrarium.saturation_specific_humidity_vapor(thermodyn_constants, air_temperature, air_density)
    relative_humidity = 0.75 # https://en.wikipedia.org/wiki/Climate_of_Norway
    specific_humidity = relative_humidity * specific_humidity_saturation
    set!(state.specific_humidity, specific_humidity)
    set!(state.ground_temperature, 13.0) # 13 °C
    set!(state.windspeed, 5.0) # 5 m/s
    compute_auxiliary!(state, model)
    @test all(isfinite.(state.skin_temperature))
    # check that skin temperature converges for a large number of iterations
    Tskin_old = deepcopy(state.skin_temperature)
    println("skin temperature at iteration 0 (default Terrarium): $(state.skin_temperature[1, 1]) net radiation: $(state.surface_net_radiation[1, 1]) Hₗ: $(state.latent_heat_flux[1, 1]) Hₛ: $(state.sensible_heat_flux[1, 1]) G: $(state.ground_heat_flux[1, 1])")
    resid = nothing
    for i in 1:10
        # compute fluxes
        Terrarium.compute_surface_energy_fluxes!(state, grid, seb, model.constants, model.atmosphere)
        # diagnose skin temperature
        Terrarium.update_skin_temperature!(state, grid, seb.skin_temperature)
        resid = maximum(abs.(state.skin_temperature - Tskin_old))
        Tskin_old = deepcopy(state.skin_temperature)
        println("skin temperature at iteration $i: $(state.skin_temperature[1, 1])  residual: $resid net radiation: $(state.surface_net_radiation[1, 1]) Hₗ: $(state.latent_heat_flux[1, 1]) Hₛ: $(state.sensible_heat_flux[1, 1]) G: $(state.ground_heat_flux[1, 1])")
    end
    @test all(state.sensible_heat_flux .> 0)
    @test all(state.ground_heat_flux .< 0)
    @test all(state.surface_net_radiation .< 0)
    @test all(resid .< sqrt(eps()))
    # Check if ground heat flux converged to gradient flux
    field_grid = get_field_grid(grid)
    Δz = Oceananigans.Δzᵃᵃᶜ(1, 1, field_grid.Nz, field_grid)
    thermal_conductivity = model.surface_energy_balance.skin_temperature.κₛ
    G_gradient = (state.ground_temperature[1, 1] - state.skin_temperature[1, 1]) / (Δz / 2 * thermal_conductivity)
    @test G_gradient ≈ state.ground_heat_flux[1, 1]

    # Cloudy and wet
    set!(state.surface_shortwave_down, 150.0)
    set!(state.surface_longwave_down, 200.0)
    air_pressure = 101_325 # standard pressure
    set!(state.air_pressure, air_pressure)
    air_temperature = 5.0 # 5 °C
    set!(state.air_temperature, air_temperature)
    air_density = Thermodynamics.air_density(thermodyn_constants, air_temperature + 273.15, air_pressure)
    specific_humidity_saturation = Terrarium.saturation_specific_humidity_vapor(thermodyn_constants, air_temperature, air_density)
    relative_humidity = 0.9 # cloudy, humid conditions
    specific_humidity = relative_humidity * specific_humidity_saturation
    set!(state.specific_humidity, specific_humidity)
    set!(state.ground_temperature, 10.0) # 10 °C
    set!(state.windspeed, 10.0) # 10 m/s
    set!(state.skin_temperature, 7.5) # initial skin temperature = ground temperature
    compute_auxiliary!(state, model)
    @test all(isfinite.(state.skin_temperature))
    # check that skin temperature converges for a large number of iterations
    Tskin_old = deepcopy(state.skin_temperature)
    for i in 1:10
        # compute fluxes
        Terrarium.compute_surface_energy_fluxes!(state, grid, seb, model.constants, model.atmosphere)
        # diagnose skin temperature
        Terrarium.update_skin_temperature!(state, grid, seb.skin_temperature)
        resid = maximum(abs.(state.skin_temperature - Tskin_old))
        Tskin_old = deepcopy(state.skin_temperature)
        println("skin temperature at iteration $i: $(state.skin_temperature[1, 1])  residual: $resid net radiation: $(state.surface_net_radiation[1, 1]) Hₗ: $(state.latent_heat_flux[1, 1]) Hₛ: $(state.sensible_heat_flux[1, 1]) G: $(state.ground_heat_flux[1, 1])")
    end
    @test all(state.sensible_heat_flux .> 0)
    @test all(state.ground_heat_flux .> 0)
    @test all(state.surface_net_radiation .> 0)
    @test all(resid .< sqrt(eps()))
end
