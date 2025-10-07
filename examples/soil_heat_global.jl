using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather

# Load land-sea mask at ~1° resolution
land_sea_frac = convert.(Float32, dropdims(Raster("data/era5-land/invariants/era5_land_land_sea_mask_N72.nc"), dims=Ti))
land_sea_frac_field = RingGrids.FullGaussianField(Matrix(land_sea_frac), input_as=Matrix)
heatmap(land_sea_frac_field)

# Load ERA-5 2 meter air temperature at ~1° resolution
Tair_raster = Raster("data/era5-land/2m_temperature/era5_land_2m_temperature_2023_N72.nc")
Tair_raster = convert.(Float32, replace_missing(Tair_raster, NaN)) .- 273.15f0
# heatmap(Tair_raster[:,:,1])

# Set up grids
land_mask = land_sea_frac_field .> 0.5 # select only grid points with > 50% land
grid = ColumnRingGrid(GPU(), ExponentialSpacing(N=30), land_mask)

# Construct input sources
forcings = InputSource(grid, rebuild(Tair_raster, name=:Tair))
# Calculate mean annual air temperature assuming hourly time resolution
Tsurf_avg = aggregate(mean, Tair_raster, (Ti(365*24), X(1), Y(1)))
heatmap(Tsurf_avg)

# Initial conditions
Tsurf_initial = Tsurf_avg[findall(land_mask)]
initializer = FieldInitializers(
    # steady-ish state initial condition calculated from initial air temperature
    temperature = (x,z) -> Tsurf_initial[Int(round(x))+1] - 0.02*z,
    # dry soil
    pore_water_ice_saturation = 1.0,
)
T_ub = PrescribedValue(:temperature, Input(:Tair, units=u"°C"))
boundary_conditions = SoilBoundaryConditions(grid, top=T_ub)
model = SoilModel(grid; initializer, boundary_conditions)
state = initialize(model, forcings)
# advance one timestep with Δt = 15 minutes
@time timestep!(state, 900.0)
# run multiple timesteps over a given time period
@time run!(state, period=Day(1), Δt=900.0)

using Oceananigans.OutputWriters: JLD2Writer, AveragedTimeInterval
using Oceananigans.Utils: days

# set up and run an Oceananigans Simulation
sim = Simulation(state; Δt=900.0, stop_iteration=100)
# sim.output_writers[:temperature] = JLD2Writer(state, get_fields(state, :temperature); filename="output/soil_model_temperature.jld2", schedule=AveragedTimeInterval(1days))
run!(sim)
