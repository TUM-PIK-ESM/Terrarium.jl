using Terrarium
using Test

import RingGrids
import Dates: Hour

@testset "run! SoilModel w/ ForwardEuler" begin
    grid = ColumnRingGrid(CPU(), Float64, ExponentialSpacing(N = 50), RingGrids.FullHEALPixGrid(16))
    model = SoilModel(grid)
    integrator = initialize(model)

    run!(integrator; steps = 2)
    @test all(isfinite.(integrator.state.temperature))

    run!(integrator; period = Hour(1))
    @test all(isfinite.(integrator.state.temperature))

    @test_throws ArgumentError run!(integrator; steps = 2, period = Hour(1))
    @test_throws ArgumentError run!(integrator)

    # test Oceananigans Simulation
    integrator = initialize(model)
    sim = Simulation(integrator; Δt = 900.0, stop_time = 3600.0)
    timestep!(sim)
    run!(sim)
    @test integrator.clock.time == 3600.0
end

@testset "run! SoilModel w/ Heun" begin
    grid = ColumnRingGrid(CPU(), Float64, ExponentialSpacing(N = 50), RingGrids.FullHEALPixGrid(16))
    model = SoilModel(grid; timestepper = Heun())
    integrator = initialize(model)

    run!(integrator; steps = 2)
    @test all(isfinite.(integrator.state.temperature))

    run!(integrator; period = Hour(1))
    @test all(isfinite.(integrator.state.temperature))

    @test_throws ArgumentError run!(integrator; steps = 2, period = Hour(1))
    @test_throws ArgumentError run!(integrator)
end
