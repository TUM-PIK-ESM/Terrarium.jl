using Terrarium
using Test

import RingGrids
import Dates: Hour
import Oceananigans: time_step!

@testset "run! SoilModel w/ ForwardEuler" begin
    grid = ColumnRingGrid(CPU(), Float64, ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16))
    model = SoilModel(; grid)
    state = initialize(model)

    run!(state; steps=2)
    @test all(isfinite.(state.state.temperature))

    run!(state; period=Hour(1))
    @test all(isfinite.(state.state.temperature))

    @test_throws ArgumentError run!(state; steps=2, period=Hour(1))
    @test_throws ArgumentError run!(state)

    # test Oceananigans Simulation
    state = initialize(model)
    sim = Simulation(state; Δt=900.0, stop_time=3600.0)
    time_step!(sim)
    run!(sim)
    @test state.clock.time == 3600.0
end 

@testset "run! SoilModel w/ Heun" begin
    grid = ColumnRingGrid(CPU(), Float64, ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16))
    model = SoilModel(; grid)
    state = initialize(model, timestepper=Heun)

    run!(state; steps=2)
    @test all(isfinite.(state.state.temperature))

    run!(state; period=Hour(1))
    @test all(isfinite.(state.state.temperature))

    @test_throws ArgumentError run!(state; steps=2, period=Hour(1))
    @test_throws ArgumentError run!(state)
end