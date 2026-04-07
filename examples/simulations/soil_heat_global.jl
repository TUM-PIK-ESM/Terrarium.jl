# # [Soil heat conduction at global scale](@id "soil_heat_global")
# Here we extend the single column soil heat conduction [example](@ref "soil_heat_column")
# to do global scale simulations, accelerated by GPU (if available).

using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets

using CairoMakie, GeoMakie

import RingGrids

input_dir = "inputs" # hide
@info "Current working directory: $(pwd())" # hide

# First we check if GPU is available and choose the architecture correspondingly.
arch = CUDA.functional() ? GPU() : CPU()
@info "Setting up simulation on $arch"

# Next, we load a land-sea mask at ~1° resolution:
land_sea_frac = convert.(Float32, dropdims(Raster(joinpath(input_dir, "era5-land_land_sea_mask_N72.nc")), dims = Ti))
land_sea_frac_field = RingGrids.FullGaussianField(Matrix(land_sea_frac), input_as = Matrix)
heatmap(land_sea_frac_field)

# Then we set up a masked [`ColumnRingGrid`](@ref), selecting only grid points
# with >50% land cover:
land_mask = land_sea_frac_field .> 0.5
grid = ColumnRingGrid(arch, Float64, ExponentialSpacing(N = 30), land_mask.grid, land_mask)

# Now we create our [`SoilModel`](@ref), this time with the default initialization:
model = SoilModel(grid)

# We impose a periodic temperature boundary condition at the surface to represent the annual cycle:
bc = PrescribedSurfaceTemperature(:T_ub, (x, t) -> 30 * sin(2π * t / (24 * 3600 * 365)))
integrator = initialize(model, ForwardEuler(), boundary_conditions = bc)

# We will do a quick check, advancing the simulation by one timestep with Δt = 15 minutes (900 seconds)
timestep!(integrator, 900.0)
@time timestep!(integrator, 900.0)

# ...then run the simulation for one day. If you have a GPU available, try changing `period = Year(1)`
# to run a full year!
@time run!(integrator, period = Day(1), Δt = 900.0)

using Oceananigans.OutputWriters: JLD2Writer, AveragedTimeInterval
using Oceananigans.Units: days

# We can also wrap the `integrator` in an Oceananigans `Simulation` which can be used to add
# callbacks and save outputs.
sim = Simulation(integrator; Δt = 900.0, stop_time = 1days)
run!(sim)
