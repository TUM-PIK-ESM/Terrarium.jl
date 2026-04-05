using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather

# run on GPU if available
arch = CUDA.functional() ? GPU() : CPU()

# Load land-sea mask at ~1° resolution
land_sea_frac = convert.(Float32, dropdims(Raster("inputs/era5-land_land_sea_mask_N72.nc"), dims = Ti))
land_sea_frac_field = RingGrids.FullGaussianField(Matrix(land_sea_frac), input_as = Matrix)
heatmap(land_sea_frac_field)

# Set up grids
land_mask = land_sea_frac_field .> 0.5 # select only grid points with > 50% land
grid = ColumnRingGrid(arch, Float64, ExponentialSpacing(N = 30), land_mask.grid, land_mask)
lon, lat = RingGrids.get_londlatds(grid.rings)

# Initial conditions
initializer = SoilInitializer(eltype(grid))
model = SoilModel(grid; initializer)
# Periodic surface temperature with annual cycle
bc = PrescribedSurfaceTemperature(:T_ub, (x, t) -> 30 * sin(2π * t / (24 * 3600 * 365)))
state = initialize(model, ForwardEuler(), boundary_conditions = bc)
# advance one timestep with Δt = 15 minutes
@time timestep!(state, 900.0)
# run multiple timesteps over a given time period
@time run!(state, period = Day(1), Δt = 900.0)

using Oceananigans.OutputWriters: JLD2Writer, AveragedTimeInterval
using Oceananigans.Units: days

# set up and run an Oceananigans Simulation
sim = Simulation(state; Δt = 900.0, stop_iteration = 100)
# sim.output_writers[:temperature] = JLD2Writer(state, get_fields(state, :temperature); filename="output/soil_model_temperature.jld2", schedule=AveragedTimeInterval(1days))
run!(sim)
