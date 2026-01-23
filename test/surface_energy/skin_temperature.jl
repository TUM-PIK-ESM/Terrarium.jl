using Terrarium
using Test

@testset "Prescribed skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    skin_temperature = PrescribedSkinTemperature()
    surface_energy_balance = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, surface_energy_balance)
    state = initialize(model)
    @test hasproperty(state.inputs, :skin_temperature)
    set!(state.skin_temperature, 1.0)
    compute_auxiliary!(state, model, surface_energy_balance)
    @test all(state.skin_temperature .≈ 1.0)
end

@testset "Implicit skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    clock = Clock(time=0.0)
    skin_temperature = ImplicitSkinTemperature()
    seb = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, seb)
    state = initialize(model)
    @test !hasproperty(state.inputs, :skin_temperature)
    @test hasproperty(state.inputs, :ground_temperature)
    set!(state.surface_shortwave_down, 300.0)
    set!(state.surface_longwave_down, 50.0)
    set!(state.specific_humidity, 0.002) # dry conditions
    set!(state.air_pressure, 101_325) # standard pressure
    set!(state.air_temperature, 2.0) # 2 °C
    set!(state.ground_temperature, 1.0) # 1 °C
    set!(state.windspeed, 1.0) # 1 m/s
    compute_auxiliary!(state, model, seb)
    @test all(isfinite.(state.skin_temperature))
    # check that skin temperature converges for a large number of iterations
    Tskin_old = deepcopy(state.skin_temperature)
    resid = nothing
    for i in 1:20
        # compute fluxes
        Terrarium.compute_surface_energy_fluxes!(state, model, seb)
        # diagnose skin temperature
        Terrarium.update_skin_temperature!(state, model, seb.skin_temperature)
        resid = maximum(abs.(state.skin_temperature - Tskin_old))
        Tskin_old = deepcopy(state.skin_temperature)
    end
    @test all(resid .< eps(Float64))
end
