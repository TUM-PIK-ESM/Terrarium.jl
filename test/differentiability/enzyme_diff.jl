using Terrarium, Enzyme
using Oceananigans: Average, Field

import SpeedyWeather.RingGrids

grid = GlobalRingGrid(CPU(), Float64, ExponentialSpacing(N=10), RingGrids.FullHEALPixGrid(16, RingGrids.Architectures.CPU()))
model = SoilModel(; grid)
sim = initialize(model)
timestep!(sim)

state = sim.state
dstate = make_zero(state)

function dostep!(state, model, ts, dt)
    timestep!(state, model, ts, dt)
    # need to use the Oceananigans reduction operators or
    # turn the Field into an Array rather than directly apply
    # reduction methods like sum, otherwise Enzyme complains.
    Tavg = Field(Average(state.temperature, dims=(1, 2, 3)))
    return Tavg[1,1,1]
end

dostep!(state, model, model.time_stepping, 1.0)

Enzyme.autodiff(set_strong_zero(set_runtime_activity(Reverse)), dostep!, Active, Duplicated(state, dstate), Const(model), Const(model.time_stepping), Const(model.time_stepping.dt))

# TODO: NaNs in the gradient, find out what's wrong 
@test_broken all(isfinite.(dstate.temperature))
