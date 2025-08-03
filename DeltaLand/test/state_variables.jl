using DeltaLand
using Test

grid = ColumnGrid(ExponentialSpacing(N=10))
model = SoilModel(; grid)
sim = initialize(model)
