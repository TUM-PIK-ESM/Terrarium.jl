using Terrarium

using CUDA
using Rasters, NCDatasets
using SpeedyWeather

using CairoMakie, GeoMakie

import RingGrids

# rg_raster = Raster(Matrix(ring_field), (X(RingGrids.get_lond(ring_grid)), Y(RingGrids.get_latd(ring_grid))))
land_sea_mask = convert.(Float32, dropdims(Raster("data/inputs/era5-land/era5_land_land_sea_mask_N145.nc"), dims=Ti))
land_sea_mask_field = RingGrids.FullGaussianField(Matrix(land_sea_mask), input_as=Matrix)
heatmap(land_sea_mask_field)

grid = ColumnRingGrid(CPU(), ExponentialSpacing(N=30), ring_grid)
# initial conditions
initializer = FieldInitializers(
    # steady-ish state initial condition for temperature
    temperature = (x,z) -> -1 - 0.02*z,
    # dry soil
    pore_water_ice_saturation = 0.0,
)
model = SoilModel(; grid, initializer)
sim = initialize(model)
@time timestep!(sim)
@time run!(sim, period=Day(10))
