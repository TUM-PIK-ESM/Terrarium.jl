using Terra, Enzyme 

import SpeedyWeather.RingGrids

grid = GlobalRingGrid(CPU(), ExponentialSpacing(N=10), RingGrids.FullHEALPixGrid(16, RingGrids.Architectures.CPU()))
initializer = FieldInitializers(temperature = (x,z) -> -1.0 - 0.01*z + exp(z/10)*sin(2π*z/10))
model = SoilModel(; grid, initializer)
sim = initialize(model)

state = sim.state
dstate = make_zero(state)

Enzyme.autodiff(set_runtime_activity(Reverse), timestep!, Const, Duplicated(state, dstate), Const(model), Const(model.time_stepping), Const(model.time_stepping.dt))
