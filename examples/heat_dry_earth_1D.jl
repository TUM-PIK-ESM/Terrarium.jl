using Terrarium

using CDSAPI
using CUDA
using Rasters

import RingGrids

ring_grid = RingGrids.FullGaussianGrid(16)

grid = ColumnRingGrid(CPU(), ExponentialSpacing(N=30), ring_grid)
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
