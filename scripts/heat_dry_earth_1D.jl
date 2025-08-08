using Terra

import SpeedyWeather.RingGrids

grid = GlobalRingGrid(GPU(), ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16, RingGrids.Architectures.GPU()))
initializer = FieldInitializers(temperature = (x,z) -> -1.0 - 0.01*z + exp(z/10)*sin(2π*z/10))
model = SoilModel(; grid, initializer)
sim = initialize(model)
@time timestep!(sim)
