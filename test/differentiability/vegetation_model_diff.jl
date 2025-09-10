using Terrarium, Enzyme
using Oceananigans: Average, Field

import RingGrids

grid = GlobalRingGrid(CPU(), Float64, ExponentialSpacing(N=10), RingGrids.FullHEALPixGrid(16, RingGrids.Architectures.CPU()))
model = VegetationModel(; grid)
sim = initialize(model)
timestep!(sim)