using Terrarium

using CUDA
using Rasters, NCDatasets

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
Tsurf_0 = Tair_raster[Ti(1)][findall(land_mask)]

# Initial conditions
initializer = FieldInitializers(
    # steady-ish state initial condition for temperature
    temperature = (x,z) -> Tsurf_0[Int(round(x))+1] - 0.02*z,
    # dry soil
    pore_water_ice_saturation = 1.0,
)
T_ub = PrescribedValue(:temperature, Input(:Tair, units=u"°C"))
boundary_conditions = SoilBoundaryConditions(grid, top=T_ub)
model = SoilModel(grid; initializer, boundary_conditions)
sim = initialize(model, forcings)
@time timestep!(sim, 900.0)
@time run!(sim, period=Day(30), dt=900.0)
