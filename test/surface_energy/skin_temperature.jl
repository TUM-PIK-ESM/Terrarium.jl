using Terrarium
using Test

@testset "Prescribed skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    skin_temperature = PrescribedSkinTemperature()
    surface_energy_balance = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, surface_energy_balance)
    model_state = initialize(model)
    state = model_state.state
    @test hasproperty(state.inputs, :skin_temperature)
    set!(state.skin_temperature, 1.0)
    compute_auxiliary!(state, model, surface_energy_balance)
    @test all(state.skin_temperature .≈ 1.0)
end

@testset "Diagnosed skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    clock = Clock(time=0.0)
    skin_temperature = ImplicitSkinTemperature()
    surface_energy_balance = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, surface_energy_balance)
    model_state = initialize(model)
    state = model_state.state
    @test !hasproperty(state.inputs, :skin_temperature)
    @test hasproperty(state.inputs, :ground_temperature)
    set!(state.surface_shortwave_down, 200.0)
    set!(state.surface_longwave_down, 100.0)
    set!(state.specific_humidity, 0.001)
    set!(state.air_pressure, 101_325)
    set!(state.air_temperature, 2.0)
    set!(state.skin_temperature, 1.0) # initial guess matching soil temperature
    set!(state.ground_temperature, 1.0)
    compute_auxiliary!(state, model, surface_energy_balance)
    @test all(state.skin_temperature - 1.96 .< 0.01)
    # check that skin temperature converges for a large number of iterations
    Tskin_old = deepcopy(state.skin_temperature)
    resid = nothing
    for i in 1:100
        compute_auxiliary!(state, model, surface_energy_balance)
        resid = maximum(abs.(state.skin_temperature - Tskin_old))
        Tskin_old = deepcopy(state.skin_temperature)
    end
    @test all(resid .< eps(Float64))
end
