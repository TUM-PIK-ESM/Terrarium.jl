using DeltaLand
using Oceananigans

import SpeedyWeather.RingGrids: FullHEALPixGrid

grid = GlobalRingGrid(ExponentialSpacing(N=50), FullHEALPixGrid(12))
model = SoilModel(; grid)
sim = initialize(model, init)
@time timestep!(sim)
