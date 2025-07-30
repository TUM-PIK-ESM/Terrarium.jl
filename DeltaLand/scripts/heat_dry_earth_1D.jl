using DeltaLand
using Oceananigans

import SpeedyWeather.RingGrids: FullHEALPixGrid

grid = GlobalRingGrid(ExponentialSpacing(N=50), FullHEALPixGrid(12))
temp_bc = FieldBoundaryConditions(grid.gridimpl, (Center,Center,Nothing); top=ValueBoundaryCondition(2.0))
init = ModelInitializer(
    :temperature => (init=(x,z) -> 1 + 0.01*abs(z), boundary_conditions=temp_bc)
)
model = SoilModel(; grid)
sim = initialize(model, init)
@time timestep!(sim)
