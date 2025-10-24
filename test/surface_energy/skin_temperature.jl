using Terrarium
using Test

@testset "Prescribed skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    skin_temperature = PrescribedSkinTemperature()
    model = SurfaceEnergyBalanceModel(grid; skin_temperature)
    model_state = initialize(model)
    state = model_state.state
    @test hasproperty(state.inputs, :skin_temperature)
    set!(state.skin_temperature, 1.0)
    compute_auxiliary!(state, model, model.radiative_fluxes)
    compute_auxiliary!(state, model, model.turbulent_fluxes)
    compute_auxiliary!(state, model, skin_temperature)
    @test all(state.skin_temperature .≈ 1.0)
end

@testset "Prognostic skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
    skin_temperature = PrognosticSkinTemperature()
    model = SurfaceEnergyBalanceModel(grid; skin_temperature)
    model_state = initialize(model)
    state = model_state.state
    @test !hasproperty(state.inputs, :skin_temperature)
    @test hasproperty(state.inputs, :ground_temperature)
    set!(state.SwIn, 200.0)
    set!(state.LwIn, 100.0)
    set!(state.specific_humidity, 0.001)
    set!(state.air_pressure, 101_325)
    set!(state.air_temperature, 2.0)
    set!(state.skin_temperature, 1.0) # initial guess matching soil temperature
    set!(state.ground_temperature, 1.0)
    compute_auxiliary!(state, model)
    Terrarium.update_skin_temperature!(state, model, skin_temperature)
    @test all(state.skin_temperature .> 1) && all(state.skin_temperature .< 2)
    # check that skin temperature converges
    Tskin_old = deepcopy(state.skin_temperature)
    resid = nothing
    for i in 1:100
        compute_auxiliary!(state, model)
        Terrarium.update_skin_temperature!(state, model, skin_temperature)
        resid = maximum(abs.(state.skin_temperature - Tskin_old))
        Tskin_old = deepcopy(state.skin_temperature)
    end
    @test all(resid .< eps(Float64))
end
