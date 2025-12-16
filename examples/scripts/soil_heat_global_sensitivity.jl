using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie
using Enzyme, Checkpointing

import RingGrids
import SpeedyWeather

# run on GPU if available
arch = CUDA.functional() ? Terrarium.GPU() : Terrarium.CPU()

# Load land-sea mask at ~1° resolution
land_sea_frac = convert.(Float32, dropdims(Raster("inputs/era5-land_land_sea_mask_N72.nc"), dims=Ti))
land_sea_frac_field = RingGrids.FullGaussianField(Matrix(land_sea_frac), input_as=Matrix)
heatmap(land_sea_frac_field)

# Set up grids
land_mask = land_sea_frac_field .> 0.5 # select only grid points with > 50% land
grid = ColumnRingGrid(arch, Float64, ExponentialSpacing(N=30), land_mask.grid, land_mask)
lon, lat = RingGrids.get_londlatds(grid.rings)

# Initial conditions
initializer = FieldInitializers(
	# steady-ish state initial condition for temperature
	temperature = (x,z) -> -1 - 0.01*z,
	# fully saturated soil pores
	saturation_water_ice = 1.0,
)
model = SoilModel(grid; initializer)
# constant surface temperature of 1°C
bcs = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator = initialize(model, ForwardEuler(), boundary_conditions=bcs)

# spin up a little
@time run!(integrator, period=Day(5), Δt=900.0)

# Enzyme prep 
scheme = Revolve(1)
dintegrator = make_zero(integrator)
N_t = 200 

# this uses checkpointing 
autodiff(set_runtime_activity(Enzyme.Reverse), run!, Const, Duplicated(integrator, dintegrator), Const(scheme), Const(N_t))

function run_sim!(integrater, N_t)
    run!(integrater, steps=N_t, Δt=900.0)
    return nothing
end

# no checkpointing 
N_t = 1 
autodiff(set_runtime_activity(Enzyme.Reverse), run_sim!, Const, Duplicated(integrator, dintegrator), Const(N_t))

    
