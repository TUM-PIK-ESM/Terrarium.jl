# test that the run! function runs the full model (on CPU)
import Terrarium.RingGrids
import Dates: Hour

@testset "run!" begin
    grid = GlobalRingGrid(CPU(), ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16, RingGrids.Architectures.CPU()))
    model = SoilModel(; grid)
    sim = initialize(model)

    run!(sim; steps=2)
    @test all(isfinite.(sim.state.temperature))

    run!(sim; period=Hour(1))
    @test all(isfinite.(sim.state.temperature))

    @test_throws ArgumentError run!(sim; steps=2, period=Hour(1))
    @test_throws ArgumentError run!(sim)
end 