using Terrarium

import SpeedyWeather.RingGrids

grid = GlobalRingGrid(GPU(), ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16, RingGrids.Architectures.GPU()))
# temperature initial condition
initializer = VarInitializer(:temperature, bounday_conditions=(top=ValueBoundaryCondition(1.0),)) do x, z
    -1 - 0.1*z + exp(z)*sin(2π*z)
end
model = SoilModel(; grid, initializer)
sim = initialize(model)
@time timestep!(sim)
