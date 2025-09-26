using Terrarium

using CUDA
using Rasters, NCDatasets

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather

# load land-sea mask
land_sea_frac = convert.(Float32, dropdims(Raster("data/inputs/era5-land/land_sea_mask/era5_land_land_sea_mask_N72.nc"), dims=Ti))
land_sea_frac_field = RingGrids.FullGaussianField(Matrix(land_sea_mask), input_as=Matrix)
heatmap(land_sea_frac_field)

# load ERA-5 2 meter air temperature
air_temperature = convert.(Float32, dropdims(Raster("data/inputs/era5-land/2m_temperature/era5_land_2m_temperature_2023_N72.nc"), dims=Ti))
air_temperature_field = RingGrids.FullGaussianField(Matrix(air_temperature), input_as=Matrix)
heatmap(air_temperature_field)

land_mask = land_sea_frac_field .> 0.5 # select only grid points with > 50% land
grid = ColumnRingGrid(GPU(), ExponentialSpacing(N=30), land_mask)
# initial conditions
initializer = FieldInitializers(
    # steady-ish state initial condition for temperature
    temperature = (x,z) -> -1 - 0.02*z,
    # dry soil
    pore_water_ice_saturation = 0.0,
)
boundary_conditions = SoilBoundaryConditions(
    top = (temperature=ValueBoundaryCondition((x, z, Tair) -> Tair, field_dependencies=(:Tair,)),)
)
model = SoilModel(; grid, initializer)
sim = initialize(model)
@time timestep!(sim)
@time run!(sim, period=Day(10))
