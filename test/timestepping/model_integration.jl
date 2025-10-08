import RingGrids
import Dates: Hour

@testset "run! SoilModel w/ ForwardEuler" begin
    grid = ColumnRingGrid(CPU(), Float64, ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16))
    model = SoilModel(; grid, time_stepping=ForwardEuler())
    sim = initialize(model)

    run!(sim; steps=2)
    @test all(isfinite.(sim.state.temperature))

    run!(sim; period=Hour(1))
    @test all(isfinite.(sim.state.temperature))

    @test_throws ArgumentError run!(sim; steps=2, period=Hour(1))
    @test_throws ArgumentError run!(sim)
end 