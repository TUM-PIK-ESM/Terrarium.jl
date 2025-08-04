using DeltaLand

import SpeedyWeather.RingGrids: FullHEALPixGrid

grid = GlobalRingGrid(ExponentialSpacing(N=50), FullHEALPixGrid(12))
initializer = FieldInitializers(temperature = (x,y,z) -> 1.0 - 0.01*z + exp(z/10)*sin(2π*z))
model = SoilModel(; grid, initializer)
sim = initialize(model)
@time timestep!(sim)
