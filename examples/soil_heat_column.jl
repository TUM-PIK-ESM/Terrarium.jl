using Terrarium

import CairoMakie as Makie

grid = ColumnGrid(ExponentialSpacing())
initializer = FieldInitializers(
    # steady-ish state initial condition for temperature
    temperature = (x,z) -> -1 - 0.01*z,
    # fully saturated soil pores
    saturation_water_ice = 1.0,
)
model = SoilModel(grid; initializer)
# constant surface temperature of 1°C
bcs = PrescribedSurfaceTemperature(:T_ub, 1.0)
driver = initialize(model, ForwardEuler, boundary_conditions=bcs)
# test one timestep
@time timestep!(driver)
# run simulation forward for a set period of time
run!(driver, period=Day(10))

T = interior(driver.state.temperature)[1,1,:]
f = interior(driver.state.liquid_water_fraction)[1,1,:]
zs = znodes(driver.state.temperature)
# Plot temperature and liquid fraction profiles
Makie.scatterlines(T, zs)
Makie.scatterlines(f, zs)