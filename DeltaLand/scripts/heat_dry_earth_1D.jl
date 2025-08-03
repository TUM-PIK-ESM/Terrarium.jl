using DeltaLand
using Oceananigans

import SpeedyWeather.RingGrids: FullHEALPixGrid

grid = GlobalRingGrid(ExponentialSpacing(N=50), FullHEALPixGrid(12))
temp_bc = FieldBoundaryConditions(grid.gridimpl, (Center,Center,Nothing); top=ValueBoundaryCondition(2.0))
model = SoilModel(; grid)
sim = initialize(model, init)
@time timestep!(sim)
