using Terrarium
using CUDA

import SpeedyWeather.RingGrids

grid = GlobalRingGrid(GPU(), ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16, RingGrids.Architectures.GPU()))
# initial conditions
initializer = Initializers(
    # steady-ish state initial condition for temperature
    VarInitializer((x,z) -> -1 - 0.01*z, :temperature),
    # dry soil
    VarInitializer(1.0, :pore_water_ice_saturation),
)
boundary_conditions = SoilBoundaryConditions(grid, top=(temperature=ValueBoundaryCondition(2.0),))
model = SoilModel(; grid, initializer, boundary_conditions)
sim = initialize(model)
@time timestep!(sim)
@time run!(sim, period=Day(10))
