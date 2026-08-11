# # [Global land simulation forced by ERA5-Land](@id land_global_era5)
#
# This example configures a fully-coupled global [`LandModel`](@ref) at ~1° resolution (N72, 72
# Gaussian rings) and drives it with one year of ERA5-Land reanalysis. The model couples
#
# * a [`SoilEnergyWaterCarbon`](@ref) column solving for transient heat conduction in all soil columns,
# * a single-layer snowpack ([`SingleLayerSnow`](@ref)),
# * prescribed vegetation ([`PrescribedVegetation`](@ref)) whose leaf area index is imposed from the
#   ERA5-Land LAI climatology, and
# * the surface energy balance and surface hydrology schemes that couple them to the atmosphere.
#
# The near-surface atmospheric state (temperature, humidity, pressure, wind, precipitation, and
# downwelling radiation) is supplied as time-varying inputs from the [`ERA5LandForcings`](@ref)
# asset via the [`InputSource`](@ref) interface, which linearly interpolates each field in time.
#
# !!! warning "Illustrative configuration"
#     This is a demonstration of the coupled-model setup, not a validated land-surface product. The
#     ERA5-Land accumulated fluxes (radiation and precipitation) are de-accumulated to hourly rates
#     as described below; the only approximation is that the small amount of flux falling exactly on
#     the nightly accumulation reset is discarded. The script is intended to be run directly and is
#     *not* executed during the documentation build, as it downloads several gigabytes of forcing
#     data and integrates the full global grid.

ENV["TERRARIUM_DEBUG"] = true

using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Rasters: Ti

using CairoMakie, GeoMakie

import RingGrids

import DisplayAs #hide

# Run on GPU if available, otherwise CPU.
arch = CUDA.functional() ? GPU() : CPU()
NF = Float32
@info "Setting up global land simulation on $arch" #hide

# ## Grid and land mask
# As in the [global soil heat example](@ref soil_heat_global), we build a masked
# [`ColumnRingGrid`](@ref) at ~1° resolution (72 Gaussian rings), keeping only grid points with more
# than 50% land cover. The land-sea mask is loaded from the ERA5-Land invariants asset and regridded
# onto the model grid. We use an exponentially-spaced soil column so that the near-surface layers,
# where the diurnal and seasonal signals are strongest, are resolved most finely. The minimum layer
# thickness is set to a relatively coarse 10 cm for now to maintain numerical stability. Future Development
# of improved timestepping methods will relax this constraint.
model_rings = on_architecture(arch, RingGrids.FullGaussianGrid(72))
land_sea_frac_native = RingGrids.Field(arch, ERA5LandInvariants(), "lsm"; NF)
land_sea_frac = RingGrids.interpolate(model_rings, land_sea_frac_native)
land_mask = land_sea_frac .> 0.5
Nz = 30 # number of soil layers
grid = ColumnRingGrid(arch, NF, ExponentialSpacing(N = Nz, Δz_min = 0.1), land_mask)

# ## Meteorological forcings from ERA5-Land
# The [`ERA5LandForcings`](@ref) asset holds one year of hourly ERA5-Land fields already regridded to
# the N72 grid. The full file is large (several GB), so we open each variable lazily and materialize
# only the time slices needed at each step. Missing (ocean) values are replaced with `NaN`.
load_forcing(name) = Terrarium.load_asset(ERA5LandForcings(), name; NF)

t2m = load_forcing("t2m")     # 2 m air temperature (K)
d2m = load_forcing("d2m")     # 2 m dewpoint temperature (K)
sp = load_forcing("sp")       # surface pressure (Pa)
u10 = load_forcing("u10")     # 10 m eastward wind component (m/s)
v10 = load_forcing("v10")     # 10 m northward wind component (m/s)
tp = load_forcing("tp")       # total precipitation (m, cumulative since 00 UTC)
sf = load_forcing("sf")       # snowfall, water equivalent (m, cumulative since 00 UTC)
ssrd = load_forcing("ssrd")   # downwelling shortwave radiation (J/m², cumulative since 00 UTC)
strd = load_forcing("strd")   # downwelling longwave radiation (J/m², cumulative since 00 UTC)

# The model consumes inputs in its own units, so we apply the necessary conversions as lazy
# broadcasts over the rasters. Only the slices touched by the time interpolation are ever evaluated,
# and because the ERA5-Land fields are already stored in `Float32` these single-level arithmetic
# broadcasts preserve the numeric type expected by the model grid.
#
# The *instantaneous* fields (temperature, dewpoint, pressure, wind) need only unit conversions.
# Near-surface specific humidity is not provided directly by ERA5-Land, so we derive it from the
# dewpoint temperature and surface pressure with [`Terrarium.dewpoint_specific_humidity`](@ref),
# using a standalone set of thermodynamic constants.
air_pressure = sp                     # Pa
wind_u = u10                          # m/s
wind_v = v10                          # m/s
constants = ThermodynamicConstants(NF)
air_temperature = Terrarium.kelvin_to_celsius.(constants, t2m)   # K -> °C
d2m_C = Terrarium.kelvin_to_celsius.(constants, d2m)
specific_humidity = Terrarium.dewpoint_specific_humidity.(constants, d2m_C, sp)

# The *accumulated* fields (radiation and precipitation) are cumulative since 00 UTC and reset each
# day, so dividing the raw values by the output interval would badly overestimate the flux. We
# instead recover the instantaneous hourly rate as the positive time difference between consecutive
# records divided by the accumulation window; `max(·, 0)` discards the daily reset, which falls at
# night where the radiative and (typically) precipitation fluxes are negligible. `deaccumulate`
# returns a lazy raster spanning the records from the second hour onward.
function deaccumulate(accumulated, accumulation_window = NF(3600))
    n = size(accumulated, Ti)
    ## the forcing rasters are dimensioned (X, Y, Ti), so the time axis is the third dimension
    current = @view accumulated[Ti(2:n)]
    previous = @view accumulated[Ti(1:(n - 1))]
    rate = max.(parent(current) .- parent(previous), zero(NF)) ./ accumulation_window
    return rebuild(previous; data = rate)
end

shortwave_down = deaccumulate(ssrd)                                           # J/m² -> W/m²
longwave_down = deaccumulate(strd)                                            # J/m² -> W/m²
total_precipitation = deaccumulate(tp)                                        # m -> m/s
snowfall = deaccumulate(sf)                                                   # m -> m/s
rainfall = max.(total_precipitation .- snowfall, zero(NF))                    # liquid = total - snow

# ## Prescribed LAI climatology
# The `ERA5LandLeafAreaIndex` asset holds a daily LAI climatology (1980–2010) at the native 0.1°
# resolution; `cycle = true` repeats it every year over the whole simulation. We use the
# high-vegetation LAI (`lai_hv`).
lai_asset = ERA5LandLeafAreaIndex()
lai_highveg = Terrarium.load_asset(lai_asset, "lai_hv"; NF, fill_value = zero(NF))

# ## Coupled land model
# We assemble the [`LandModel`](@ref) from its process components. The soil simulates transient heat conduction
# with static soil hydrology (for now); the vegetation is prescribed (LAI imposed, no prognostic carbon pool);
# a single-layer snowpack accumulates snowfall and couples to the surface energy balance. The `PrescribedAtmosphere`
# uses the [`WindVelocity`](@ref) parameterization so that the ERA5 `u`/`v` wind components can be mapped directly to
# the `windspeed` input required for the surface energy balance.
soil = SoilEnergyWaterCarbon(NF)
vegetation = PrescribedVegetation(NF)
snow = SingleLayerSnow(NF)
atmosphere = PrescribedAtmosphere(NF; wind = WindVelocity())
land = LandModel(grid; soil, vegetation, snow, atmosphere)
@show variables(land)

# ## Assembling the input sources
# Each converted field becomes a time-varying [`InputSource`](@ref) bound to the model grid. The
# meteorological forcings already live on the model's N72 rings, so no spatial regridding is needed;
# the LAI climatology is regridded from its native grid. We anchor every source to a common `reftime`
# (the first record of the meteorological file) so the de-accumulated fluxes — whose time axis starts
# one hour later — stay aligned in time with the instantaneous fields; simulation time `t = 0` then
# corresponds to the start of the forcing year.
reftime = first(dims(t2m, Ti))
inputs = InputSources(
    InputSource(grid, air_temperature; name = :air_temperature, units = u"°C", reftime),
    InputSource(grid, air_pressure; name = :air_pressure, units = u"Pa", reftime),
    InputSource(grid, wind_u; name = :wind_u, units = u"m/s", reftime),
    InputSource(grid, wind_v; name = :wind_v, units = u"m/s", reftime),
    InputSource(grid, specific_humidity; name = :specific_humidity, units = u"kg/kg", reftime),
    InputSource(grid, rainfall; name = :rainfall, units = u"m/s", reftime),
    InputSource(grid, snowfall; name = :snowfall, units = u"m/s", reftime),
    InputSource(grid, shortwave_down; name = :surface_shortwave_down, units = u"W/m^2", reftime),
    InputSource(grid, longwave_down; name = :surface_longwave_down, units = u"W/m^2", reftime),
    InputSource(grid, lai_highveg; source_grid = Terrarium.native_grid(lai_asset), name = :leaf_area_index, cycle = true),
)

# ## Initial conditions
# The soil moisture starts from a variably-saturated profile that approaches saturation with depth,
# and the surface starts snow-free so the winter forcing builds up the snowpack. The soil temperature
# is anchored to the local ERA5 near-surface air temperature: we initialize the model once (which
# loads the forcings at `t = 0`), read the resulting per-column air temperature from the state, and
# re-initialize with a soil temperature profile built from it plus a mild geothermal gradient. Reading
# the air temperature back from the state guarantees it is ordered consistently with the land columns.
initializers = (
    saturation_water_ice = (x, z) -> min(one(NF), NF(0.6) - NF(0.05) * z),
)
integrator = initialize(land; inputs, initializers)

# Helper to plot the top (surface) layer of a 3D field on the ring grid.
plot_surface(field; kwargs...) = heatmap(RingGrids.Field(field, grid)[:, end]; kwargs...)

# The initial surface soil temperature and prescribed (mid-January) leaf area index:
fig = plot_surface(integrator.state.temperature, title = "Initial surface soil temperature (°C)")
DisplayAs.PNG(fig) #hide

fig = plot_surface(integrator.state.leaf_area_index, title = "Leaf area index (Jan)", colorrange = (0, 6))
DisplayAs.PNG(fig) #hide

# ## Run through the first three months
# We will advance the coupled model for one day. Note that, due to  both timestepping restrictions and I/O overhead,
# the simulation is currently very slow. This will improve in the near future!
Terrarium.initialize!(integrator)
@profview @time timestep!(integrator)
run!(integrator, period = Day(1), Δt = Minute(15), show_progress = true)

# Let's look at the results:
# First we'll inspect the topsoil temperature after one month of forcing.
fig = plot_surface(integrator.state.temperature, title = "Surface soil temperature")
DisplayAs.PNG(fig) #hide

# We can also see that snow has built up over the high latitudes and elevated terrain of the winter (northern) hemisphere.
fig = plot_surface(integrator.state.snow_water_equivalent, title = "Snow water equivalent (m)")
DisplayAs.PNG(fig) #hide

# The surface latent heat flux reflects the coupled control of soil moisture, vegetation, and the
# atmospheric state on evapotranspiration.
fig = plot_surface(integrator.state.latent_heat_flux, title = "Latent heat flux (W/m²)")
DisplayAs.PNG(fig) #hide
