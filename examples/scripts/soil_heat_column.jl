using Terrarium

import CairoMakie as Makie

# run on GPU if available
arch = CUDA.functional() ? GPU() : CPU()
# create single column grid
grid = ColumnGrid(arch, Float32, ExponentialSpacing(N = 10))
# initializer for soil model with quasi thermal steady-state
initializer = SoilInitializer(
    energy = QuasiThermalSteadyState(T₀ = -1.0)
)
model = SoilModel(grid; initializer)
# constant surface temperature of 1°C
boundary_conditions = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator = initialize(model, ForwardEuler(); boundary_conditions)
# test one timestep
@time timestep!(integrator)
# run simulation forward for a set period of time
@time run!(integrator, period = Day(3))

T = interior(integrator.state.temperature)[1, 1, :]
f = interior(integrator.state.liquid_water_fraction)[1, 1, :]
zs = znodes(integrator.state.temperature)
# Plot temperature and liquid fraction profiles
Makie.scatterlines(T, zs)
Makie.scatterlines(f, zs)
