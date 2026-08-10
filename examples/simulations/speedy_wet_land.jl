using Terrarium

using CUDA
using Dates
using NumericalEarth
using NumericalEarth.DataWrangling
using NumericalEarth.DataWrangling.SoilGrids
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather as Speedy

# Choose architecture based on available hardware
arch = CUDA.functional() ? GPU() : CPU()

# Build the SpeedyWeather ring + spectral grid
speedy_arch = RingGrids.Architectures.architecture(arch)
spectral_grid = Speedy.SpectralGrid(truncation = 64, architecture = speedy_arch)
ring_grid = spectral_grid.grid

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

function Terrarium.InputSources(dataset::SoilGrids2, grid::ColumnRingGrid, horizons = (Symbol(:horizon, i) for i in 1:6); name = nameof(typeof(dataset)), verbose = true)
    soilgrids_vars = (:sand_fraction, :silt_fraction, :clay_fraction, :bulk_density)
    metadataset = MetadataSet(soilgrids_vars...; dataset)
    arch = RingGrids.architecture(grid.rings)
    soilgrids_inputs = []
    for (idx, horizon) in enumerate(horizons)
        layer_inputs = Dict()
        for var in soilgrids_vars
            verbose && @info "Loading input data for $var on $horizon"
            var_field = Field(getproperty(metadataset, var))
            ring_field = RingGrids.on_architecture(arch, RingGrids.FullClenshawField(interior(var_field)[:, (end - 1):-1:2, end - idx + 1], input_as = Matrix))
            target_field = RingGrids.Field(grid.rings)
            RingGrids.interpolate!(target_field, ring_field)
            layer_inputs[var] = InputSource(grid, Field(target_field, grid); name = horizon => var)
        end
        # Ensure that mineral texture components with each horizon sum to unity
        Terrarium.normalize_texture!(layer_inputs[:sand_fraction].field, layer_inputs[:silt_fraction].field, layer_inputs[:clay_fraction].field)
        append!(soilgrids_inputs, values(layer_inputs))
    end
    return InputSources(name, soilgrids_inputs...)
end

# Build the Terrarium land model on the matching ring grid with appropriate land-sea mask
Nz = 30
Δz_min = 0.05
land_sea_mask = Speedy.EarthLandSeaMask(spectral_grid)
Speedy.load_mask!(land_sea_mask)
Makie.heatmap(on_architecture(CPU(), land_sea_mask.land_fraction))
land_grid = ColumnRingGrid(arch, Float32, ExponentialSpacing(; N = Nz, Δz_min), land_sea_mask.land_fraction .> 0)
porosity = SoilPorositySURFEX(eltype(land_grid))
strat = SoilGridsStratigraphy(eltype(land_grid); porosity)
soil = SoilEnergyWaterCarbon(eltype(land_grid); strat)
terrarium_model = Terrarium.LandModel(land_grid; vegetation = nothing, soil)

# Wrap the Terrarium model as a SpeedyWeather land component
soilgrids_inputs = InputSources(SoilGrids2(), land_grid)
land = Speedy.LandModel(
    spectral_grid,
    terrarium_model;
    inputs = soilgrids_inputs,
    initializers = (temperature = initial_soil_temperature(land_grid),),
    Δt = Minute(5)
)

# Build the coupled PrimitiveWetModel
surface_heat_flux = Speedy.SurfaceHeatFlux(land.spectral_grid, land = Speedy.PrescribedLandHeatFlux())
surface_humidity_flux = Speedy.SurfaceHumidityFlux(land.spectral_grid, land = Speedy.PrescribedLandHumidityFlux())
output = Speedy.NetCDFOutput(land.spectral_grid, Speedy.PrimitiveDryModel; interval = Hour(3), nlayers_soil = land.geometry.nlayers, path = "outputs/")
time_stepping = Speedy.Leapfrog(land.spectral_grid, Δt_at_T32 = Minute(15))
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
sim_coupled = @time Speedy.initialize!(primitive_wet_coupled)
period = Month(6)
@info "Running simulation for $period"
@time Speedy.run!(sim_coupled, period = period, output = true)

# The Terrarium state lives inside SpeedyWeather's Variables tree
land_state = sim_coupled.variables.prognostic.land.terrarium
Terrarium.checkfinite!(land_state.prognostic)

# Here we declare some plotting helpers that ensure all Fields/arrays are transferred to CPU before attempting to plot:
plot_land_field(field, z_idx=1; kwargs...) = heatmap(RingGrids.Field(CPU(), interior(field)[:, :, z_idx], land_grid); kwargs...)
plot_speedy_field(field; kwargs...) = heatmap(on_architecture(CPU(), field); kwargs...)

# Land variables (use the SpeedyWeather-owned Terrarium state)
Tsoil_fig = plot_land_field(land_state.temperature, Nz - 1, title = "Soil temperature (10 cm depth)", size = (800, 400))
Tsurf_fig = plot_land_field(land_state.skin_temperature, title = "Skin temperature", size = (800, 400))
sat_fig = plot_land_field(land_state.saturation_water_ice, title = "Soil saturation water + ice", size = (800, 400))
Hs_fig = plot_land_field(land_state.sensible_heat_flux, title = "Sensible heat flux", size = (800, 400))
Hl_fig = plot_land_field(land_state.latent_heat_flux, title = "Latent heat flux", size = (800, 400))
E_fig = plot_land_field(land_state.evaporation_ground, title = "Evaporation (bare ground)", size = (800, 400))
swe_fig = plot_land_field(land_state.snow_water_equivalent, title = "Snow water equivalent", size = (800, 400))
scf_fig = plot_land_field(land_state.snow_cover_fraction, title = "Snow cover fraction", size = (800, 400))
snowT_fig = plot_land_field(land_state.snow_temperature, title = "Snow temperature", size = (800, 400))

# Atmosphere variables
Tair_fig = plot_speedy_field(sim_coupled.variables.parameterizations.surface_air_temperature .- 273.15, title = "Air temperature", size = (800, 400))
prcp_fig = plot_speedy_field(sim_coupled.variables.parameterizations.rain_rate, title = "Rainfall", size = (800, 400))
snow_fig = plot_speedy_field(sim_coupled.variables.parameterizations.snow_rate, title = "Snowfall", size = (800, 400))
pres_fig = plot_speedy_field(exp.(sim_coupled.variables.grid.pressure), title = "Surface pressure", size = (800, 400))
srad_fig = plot_speedy_field(sim_coupled.variables.parameterizations.surface_shortwave_down, title = "Surface shortwave down", size = (800, 400))

# Pick a point somewhere in the mid-latitudes and plot the vertical soil layers
T = on_architecture(CPU(), interior(land_state.temperature)[2000, 1, :])
sat = on_architecture(CPU(), interior(land_state.saturation_water_ice)[2000, 1, :])
f = on_architecture(CPU(), interior(land_state.liquid_water_fraction)[2000, 1, :])
zs = on_architecture(CPU(), znodes(land_state.temperature))
Makie.scatterlines(T, zs)
Makie.scatterlines(sat, zs)
Makie.scatterlines(f, zs)

# Make some animations!
sim_copuled_cpu = on_architecture(CPU(), sim_coupled)
Speedy.animate(sim_copuled_cpu, output_file = "plots/speedy_terrarium_tair_animation.mp4", variable = "temp", framerate = 12, level = spectral_grid.nlayers, transient_timesteps = 1)
Speedy.animate(sim_copuled_cpu, output_file = "plots/speedy_terrarium_tsoil_animation.mp4", variable = "st", framerate = 12, level = 1, transient_timesteps = 1)
Speedy.animate(sim_copuled_cpu, output_file = "plots/speedy_terrarium_shf_animation.mp4", variable = "shf", framerate = 12, transient_timesteps = 1)
Speedy.animate(sim_copuled_cpu, output_file = "plots/speedy_terrarium_shuf_animation.mp4", variable = "shuf", framerate = 12, transient_timesteps = 1)
Speedy.animate(sim_copuled_cpu, output_file = "plots/speedy_terrarium_sd_animation.mp4", variable = "sd", framerate = 12, transient_timesteps = 1)
