using Terrarium
using Test

import RingGrids
import Dates: Hour
import Oceananigans: time_step!

@testset "run! SoilModel w/ ForwardEuler" begin
    grid = ColumnRingGrid(CPU(), Float64, ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16))
    model = SoilModel(grid)
    driver = initialize(model, ForwardEuler())

    run!(driver; steps=2)
    @test all(isfinite.(driver.state.temperature))

    run!(driver; period=Hour(1))
    @test all(isfinite.(driver.state.temperature))

    @test_throws ArgumentError run!(driver; steps=2, period=Hour(1))
    @test_throws ArgumentError run!(driver)

    # test Oceananigans Simulation
    driver = initialize(model, ForwardEuler())
    sim = Simulation(driver; Δt=900.0, stop_time=3600.0)
    time_step!(sim)
    run!(sim)
    @test driver.clock.time == 3600.0
end 

@testset "run! SoilModel w/ Heun" begin
    grid = ColumnRingGrid(CPU(), Float64, ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16))
    model = SoilModel(grid)
    driver = initialize(model, Heun())

    run!(driver; steps=2)
    @test all(isfinite.(driver.state.temperature))

    run!(driver; period=Hour(1))
    @test all(isfinite.(driver.state.temperature))

    @test_throws ArgumentError run!(driver; steps=2, period=Hour(1))
    @test_throws ArgumentError run!(driver)
end