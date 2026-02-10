using Terrarium
using Test

@testset "Prescribed skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))
    skin_temperature = PrescribedSkinTemperature(eltype(grid))
    seb = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, seb)
    state = initialize(model)
    @test hasproperty(state.inputs, :skin_temperature)
    set!(state.skin_temperature, 1.0)
    compute_auxiliary!(state, grid, skin_temperature)
    @test all(state.skin_temperature .≈ 1.0)
end

@testset "Implicit skin temperature" begin
    grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))
    clock = Clock(time = 0.0)
    skin_temperature = ImplicitSkinTemperature()
    seb = SurfaceEnergyBalance(Float64; skin_temperature)
    model = SurfaceEnergyModel(grid, seb)
    state = initialize(model)
    @test !hasproperty(state.inputs, :skin_temperature)
    @test hasproperty(state.inputs, :ground_temperature)
    set!(state.surface_shortwave_down, 300.0) # sunny conditions
    set!(state.surface_longwave_down, 50.0)
    set!(state.specific_humidity, 0.002) # dry conditions
    set!(state.air_pressure, 101_325) # standard pressure
    set!(state.air_temperature, 10.0) # 10 °C
    set!(state.ground_temperature, 2.0) # 2 °C
    set!(state.windspeed, 1.0) # 1 m/s
    compute_auxiliary!(state, grid, skin_temperature, seb)
    @test all(isfinite.(state.skin_temperature))
    # check that skin temperature converges for a large number of iterations
    Tskin_old = deepcopy(state.skin_temperature)
    resid = nothing
    for i in 1:5
        # compute fluxes
        Terrarium.compute_surface_energy_fluxes!(state, grid, seb, model.atmosphere, model.constants)
        # diagnose skin temperature
        Terrarium.update_skin_temperature!(state, grid, seb.skin_temperature)
        resid = maximum(abs.(state.skin_temperature - Tskin_old))
        Tskin_old = deepcopy(state.skin_temperature)
        println("skin temperature residual at iteration $i: $resid")
    end
    @test all(resid .< sqrt(eps()))
end
