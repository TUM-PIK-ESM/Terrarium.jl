using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather as Speedy

# Choose architecture based on available hardware
arch = CPU() #CUDA.functional() ? GPU() : CPU()

# Build the SpeedyWeather ring + spectral grid
speedy_arch = RingGrids.Architectures.architecture(arch)
spectral_grid = Speedy.SpectralGrid(trunc = 31, architecture = speedy_arch)
ring_grid = spectral_grid.grid

# Build the Terrarium land model on the matching ring grid
Nz = 30
Δz_min = 0.05
land_grid = ColumnRingGrid(arch, Float32, ExponentialSpacing(; N = Nz, Δz_min), ring_grid)
soil = SoilEnergyWaterCarbon(eltype(land_grid), hydrology = SoilHydrology(eltype(land_grid)))
terrarium_model = Terrarium.LandModel(land_grid; vegetation = nothing, soil)

# To make the simulation a bit more interesting, we will use spatially periodic initial and boundary conditions.
# The climatology will be determined by latitude with a maximum of 20 °C at the equator and minimum of -20°C at
# the poles.
mean_annual_temperature(lat) = 20 - abs(40 * sin(lat)) # maximum at equator

function initial_soil_temperature(grid)
    arch = Terrarium.architecture(grid)
    _, grid_lat = RingGrids.get_lonlats(grid.rings) # in radians
    function init(x, z)
        latᵢ = grid_lat[round(Int, x)]
        T₀ = mean_annual_temperature(latᵢ)
        T = T₀ - 0.05 * z
        return T
    end
    return init
end

# Wrap the Terrarium model as a SpeedyWeather land component
land = Speedy.LandModel(
    spectral_grid,
    terrarium_model;
    initializers = (temperature = initial_soil_temperature(land_grid),),
    Δt = Minute(5)
)

# Build the coupled PrimitiveWetModel
land_sea_mask = Speedy.RockyPlanetMask(land.spectral_grid)
surface_heat_flux = Speedy.SurfaceHeatFlux(land.spectral_grid, land = Speedy.PrescribedLandHeatFlux())
surface_humidity_flux = Speedy.SurfaceHumidityFlux(land.spectral_grid, land = Speedy.PrescribedLandHumidityFlux())
output = Speedy.NetCDFOutput(land.spectral_grid, Speedy.PrimitiveDryModel; interval = Hour(3), nlayers_soil = land.geometry.nlayers, path = "outputs/")
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
Speedy.add!(primitive_wet_coupled.output, Speedy.SnowDepthOutput())

# Run the coupled simulation
@info "Initializing coupled simulation"
sim_coupled = Speedy.initialize!(primitive_wet_coupled)
period = Month(6)
@info "Running simulation for $period"
@time Speedy.run!(sim_coupled, period = period, output = true)

# The Terrarium state lives inside SpeedyWeather's Variables tree
land_state = sim_coupled.variables.prognostic.land.terrarium
Terrarium.checkfinite!(land_state.prognostic)

# Land variables (use the SpeedyWeather-owned Terrarium state)
Tsoil_fig = heatmap(RingGrids.Field(interior(land_state.temperature)[:, 1, end - 3], grid), title = "", size = (800, 400))
Tsurf_fig = heatmap(RingGrids.Field(interior(land_state.skin_temperature)[:, 1], grid), title = "", size = (800, 400))
Hs_fig = heatmap(RingGrids.Field(interior(land_state.sensible_heat_flux)[:, 1], grid), title = "", size = (800, 400))
Hl_fig = heatmap(RingGrids.Field(interior(land_state.latent_heat_flux)[:, 1], grid), title = "", size = (800, 400))
E_fig = heatmap(RingGrids.Field(interior(land_state.evaporation_ground)[:, 1], grid), title = "", size = (800, 400))
swe_fig = heatmap(RingGrids.Field(interior(land_state.snow_water_equivalent)[:, 1, end], grid), title = "", size = (800, 400))
scf_fig = heatmap(RingGrids.Field(interior(land_state.snow_cover_fraction)[:, 1, end], grid), title = "", size = (800, 400))
snowT_fig = heatmap(RingGrids.Field(interior(land_state.snow_temperature)[:, 1, end], grid), title = "", size = (800, 400))
sat_fig = heatmap(RingGrids.Field(interior(land_state.saturation_water_ice)[:, 1, end], grid), title = "", size = (800, 400))

# Atmosphere variables
Tair_fig = heatmap(sim_coupled.variables.parameterizations.surface_air_temperature .- 273.15, title = "Air temperature", size = (800, 400))
prcp_fig = heatmap(sim_coupled.variables.parameterizations.rain_rate, title = "Rainfall", size = (800, 400))
snow_fig = heatmap(sim_coupled.variables.parameterizations.snow_rate, title = "Snowfall", size = (800, 400))
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

# Make some animations!
Speedy.animate("outputs/run_0001/output.nc", output_file = "plots/speedy_terrarium_tair_animation.mp4", variable = "temp", framerate = 12, level = spectral_grid.nlayers, transient_timesteps = 2)
Speedy.animate("outputs/run_0001/output.nc", output_file = "plots/speedy_terrarium_tsoil_animation.mp4", variable = "st", framerate = 12, level = 1, transient_timesteps = 2)
Speedy.animate("outputs/run_0001/output.nc", output_file = "plots/speedy_terrarium_shf_animation.mp4", variable = "shf", framerate = 12, transient_timesteps = 2)
Speedy.animate("outputs/run_0001/output.nc", output_file = "plots/speedy_terrarium_shuf_animation.mp4", variable = "shuf", framerate = 12, transient_timesteps = 2)
Speedy.animate("outputs/run_0001/output.nc", output_file = "plots/speedy_terrarium_sd_animation.mp4", variable = "sd", framerate = 12, transient_timesteps = 2)
