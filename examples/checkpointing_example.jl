using Terrarium, Checkpointing, Enzyme

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
integrator = initialize(model, ForwardEuler(), boundary_conditions=bcs)
scheme = Revolve(1)
N_t = 100

dintegrator = make_zero(integrator)
# set a one hot seed for a sensitivity analysis of T for now 
interior(dintegrator.state.temperature)[1,1,2] = 1.0

autodiff(set_runtime_activity(Reverse), run!, Const, Duplicated(integrator, dintegrator), Const(scheme), Const(N_t))

dU = interior(dintegrator.state.internal_energy)[1,1,:]
dT = interior(dintegrator.state.temperature)[1,1,:]
zs = znodes(integrator.state.temperature)

Makie.scatterlines(dU, zs)
Makie.scatterlines(dT, zs)
