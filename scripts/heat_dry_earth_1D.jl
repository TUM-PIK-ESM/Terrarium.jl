using Terrarium
using CUDA

import SpeedyWeather.RingGrids

ring_grid = RingGrids.FullHEALPixGrid(16, RingGrids.Architectures.GPU())
grid = GlobalRingGrid(GPU(), ExponentialSpacing(N=50), ring_grid)
# initial conditions
initializer = FieldInitializers(
    # steady-ish state initial condition for temperature
    temperature = (x,z) -> -1 - 0.02*z,
    # dry soil
    pore_water_ice_saturation = 0.0,
)
model = SoilModel(; grid, initializer)
sim = initialize(model)
@time timestep!(sim)
@time run!(sim, period=Day(10))
