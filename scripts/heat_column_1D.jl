using Terrarium
using Dates

import CairoMakie as Makie

grid = ColumnGrid(ExponentialSpacing(Δz_min=0.05, Δz_max=10.0, N=30))
initializer = Initializers(
    # steady-ish state initial condition for temperature
    VarInitializer((x,z) -> -1 - 0.05*z, :temperature),
    # fully saturated soil pores
    VarInitializer(1.0, :pore_water_ice_saturation),
)
# temperature boundary condition
boundary_conditions = VarBoundaryConditions(:temperature, top=ValueBoundaryCondition(1.0))
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