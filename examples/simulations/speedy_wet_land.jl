using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather as Speedy

# Choose architecture based on available hardware
arch = CUDA.functional() ? GPU() : CPU()

# Build the SpeedyWeather ring + spectral grid
ring_grid = RingGrids.FullGaussianGrid(24)
spectral_grid = Speedy.SpectralGrid(ring_grid)

# Build the Terrarium land model on the matching ring grid
Nz = 30
Δz_min = 0.05
grid = ColumnRingGrid(CPU(), Float32, ExponentialSpacing(; N = Nz, Δz_min), ring_grid)
soil_initializer = SoilInitializer(eltype(grid))
soil = SoilEnergyWaterCarbon(eltype(grid), hydrology = SoilHydrology(eltype(grid)))
terrarium_model = Terrarium.LandModel(grid; initializer = soil_initializer, vegetation = nothing, soil)

# Wrap the Terrarium model as a SpeedyWeather land component
land = Speedy.LandModel(
    spectral_grid, terrarium_model;
    timestepper = ForwardEuler(eltype(grid)),
    Δt = 300.0,
)

# Build the coupled PrimitiveWetModel
land_sea_mask = Speedy.RockyPlanetMask(land.spectral_grid)
surface_heat_flux = Speedy.SurfaceHeatFlux(land.spectral_grid, land = Speedy.PrescribedLandHeatFlux())
surface_humidity_flux = Speedy.SurfaceHumidityFlux(land.spectral_grid, land = Speedy.PrescribedLandHumidityFlux())
output = Speedy.NetCDFOutput(land.spectral_grid, Speedy.PrimitiveDryModel; nlayers_soil = land.geometry.nlayers, path = "outputs/")
time_stepping = Speedy.Leapfrog(land.spectral_grid, Δt_at_T31 = Minute(15))
primitive_wet_coupled = Speedy.PrimitiveWetModel(
    land.spectral_grid;
    land,
    surface_heat_flux,
    surface_humidity_flux,
    land_sea_mask,
    time_stepping,
    output,
)
Speedy.add!(primitive_wet_coupled.output, Speedy.SoilTemperatureOutput())

# Run the coupled simulation
sim_coupled = Speedy.initialize!(primitive_wet_coupled)
period = Day(1)
@info "Running simulation for $period"
@time Speedy.run!(sim_coupled, period = period)

# The Terrarium state lives inside SpeedyWeather's Variables tree
land_state = sim_coupled.variables.prognostic.land.terrarium
Terrarium.checkfinite!(land_state.prognostic)

# Land variables (use the SpeedyWeather-owned Terrarium state)
Tsoil_fig = heatmap(RingGrids.Field(interior(land_state.temperature)[:, 1, end - 2], grid), title = "", size = (800, 400))
Tsurf_fig = heatmap(RingGrids.Field(interior(land_state.skin_temperature)[:, 1], grid), title = "", size = (800, 400))
Hs_fig = heatmap(RingGrids.Field(interior(land_state.sensible_heat_flux)[:, 1], grid), title = "", size = (800, 400))
Hl_fig = heatmap(RingGrids.Field(interior(land_state.latent_heat_flux)[:, 1], grid), title = "", size = (800, 400))
E_fig = heatmap(RingGrids.Field(interior(land_state.evaporation_ground)[:, 1], grid), title = "", size = (800, 400))
sat_fig = heatmap(RingGrids.Field(interior(land_state.saturation_water_ice)[:, 1, end], grid), title = "", size = (800, 400))

# Atmosphere variables
Tair_fig = heatmap(sim_coupled.variables.grid.temperature[:, 8] .- 273.15, title = "Air temperature", size = (800, 400))
pres_fig = heatmap(exp.(sim_coupled.variables.grid.pressure), title = "Surface pressure", size = (800, 400))
srad_fig = heatmap(sim_coupled.variables.parameterizations.surface_shortwave_down, title = "Surface shortwave down", size = (800, 400))

# Pick a point somewhere in the mid-latitudes and plot the soil column
T = interior(land_state.temperature)[2000, 1, :]
sat = interior(land_state.saturation_water_ice)[2000, 1, :]
f = interior(land_state.liquid_water_fraction)[2000, 1, :]
zs = znodes(land_state.temperature)

Makie.scatterlines(T, zs)
Makie.scatterlines(sat, zs)
Makie.scatterlines(f, zs)
