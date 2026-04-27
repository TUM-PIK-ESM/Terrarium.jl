using Terrarium

import RingGrids

arch_oceananigans = CPU()
arch_speedy = RingGrids.CPU()
NF = Float32
ring_grid = RingGrids.FullGaussianGrid(72, arch_speedy)
grid = ColumnRingGrid(arch_oceananigans, NF, ExponentialSpacing(N = 30), ring_grid)
model = SoilModel(grid)
bc = PrescribedSurfaceTemperature(:T_ub, NF(10.0))
inits = (temperature = NF(15.0),)
integrator = initialize(model, ForwardEuler(NF), boundary_conditions = bc, initializers = inits)
timestep!(integrator)
@time timestep!(integrator)
