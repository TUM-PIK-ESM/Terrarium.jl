using Terrarium

import SpeedyWeather.RingGrids

grid = GlobalRingGrid(GPU(), ExponentialSpacing(N=50), RingGrids.FullHEALPixGrid(16, RingGrids.Architectures.GPU()))
# initial conditions
initializer = Initializers(
    # steady-ish state initial condition for temperature
    VarInitializer((x,z) -> -1 - 0.05*z, :temperature),
    # dry soil
    VarInitializer(0.0, :pore_water_ice_saturation),
)
model = SoilModel(; grid, initializer)
sim = initialize(model)
@time timestep!(sim)
