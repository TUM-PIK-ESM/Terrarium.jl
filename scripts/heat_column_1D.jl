using Terrarium
using Dates

import CairoMakie as Makie

grid = ColumnGrid(ExponentialSpacing())
initializer = Initializers(
    # steady-ish state initial condition for temperature
    VarInitializer((x,z) -> -1 - 0.05*z, :temperature),
    # fully saturated soil pores
    VarInitializer(1.0, :pore_water_ice_saturation),
)
# temperature boundary condition
boundary_conditions = SoilBoundaryConditions(
    grid,
    top = (temperature = ValueBoundaryCondition(1.0),)
)
model = SoilModel(; grid, initializer, boundary_conditions)
sim = initialize(model)
# test one timestep
@time timestep!(sim)
# run simulation forward for a set period of time
run!(sim, period=Day(10))

T = interior(sim.state.temperature)
f = interior(sim.state.liquid_water_fraction)
zs = znodes(sim.state.temperature)
# Plot temperature and liquid fraction profiles
Makie.scatterlines(T[1,1,:], zs)
Makie.scatterlines(f[1,1,:], zs)