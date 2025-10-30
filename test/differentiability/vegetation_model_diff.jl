using Terrarium
using Enzyme

grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N=10))
model = VegetationModel(; grid)
sim = initialize(model)
timestep!(sim)

# TODO
