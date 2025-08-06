using DeltaLand
using Dates

import Oceananigans: Center
import Oceananigans.BoundaryConditions: ValueBoundaryCondition, FieldBoundaryConditions
import SpeedyWeather.RingGrids

grid = ColumnGrid(ExponentialSpacing(Δz_min=0.05, Δz_max=10.0, N=30))
initializer = FieldInitializers(temperature = (x,z) -> -1.0 - 0.1*z)
boundary_conditions = DeltaLand.FieldBoundaryConditions(temperature=FieldBoundaryConditions(grid.grid, (Center,Center,Nothing), top=ValueBoundaryCondition(0.0)))
model = SoilModel(; grid, initializer, boundary_conditions)
sim = initialize(model)
run!(sim, period=Day(30))

sim.state.internal_energy_tendency[1,1,:]

import Adapt: adapt
import CairoMakie as Makie

T = adapt(Array, sim.state.temperature)[1,1,1:end-3]
z_centers = grid.grid.z.cᵃᵃᶜ[1:end-3]
Makie.lines(z_centers, T)
