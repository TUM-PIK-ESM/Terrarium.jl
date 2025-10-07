using Terrarium

import CairoMakie as Makie

grid = ColumnGrid(ExponentialSpacing())
initializer = FieldInitializers(
    # steady-ish state initial condition for temperature
    temperature = (x,z) -> -2 - 0.01*z,
    # fully saturated soil pores
    pore_water_ice_saturation = 1.0,
)
# temperature boundary condition
boundary_conditions = SoilBoundaryConditions(
    grid,
    # Set a constant temperature of 1°C at the upper boundary
    top = (temperature = ValueBoundaryCondition(1.0),)
)
model = SoilModel(; grid, initializer, boundary_conditions)
sim = initialize(model)
# test one timestep
@time timestep!(sim)
# run simulation forward for a set period of time
run!(sim, period=Day(10))

T = interior(sim.state.temperature)[1,1,:]
f = interior(sim.state.liquid_water_fraction)[1,1,:]
zs = znodes(sim.state.temperature)
# Plot temperature and liquid fraction profiles
Makie.scatterlines(T, zs)
Makie.scatterlines(f, zs)