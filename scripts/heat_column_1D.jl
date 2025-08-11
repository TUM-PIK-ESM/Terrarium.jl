using Terrarium
using Dates

import CairoMakie as Makie

import Oceananigans.BoundaryConditions: ValueBoundaryCondition

grid = ColumnGrid(ExponentialSpacing(Δz_min=0.05, Δz_max=100.0, N=50))
initializer = FieldInitializers(
    # steady state temperature profile
    temperature = (x,z) -> -1.0 - 0.01*z,
    # dry soil
    pore_water_ice_saturation = (x,z) -> 0.0,
)
boundary_conditions = FieldBoundaryConditions(temperature=(top=ValueBoundaryCondition(0.0),))
model = SoilModel(; grid, initializer, boundary_conditions)
sim = initialize(model)
timestep!(sim)
run!(sim, period=Day(30))

Terrarium.invclosure!(sim.state, model, Terrarium.TemperatureEnergyClosure())

# TODO: Figure out how to retrieve field data without halo regions...
T = adapt(Array, sim.state.temperature)[1,1,1:end-3]
z_centers = grid.grid.z.cᵃᵃᶜ[1:end-3]
Makie.scatterlines(T[10:end], z_centers[10:end])
