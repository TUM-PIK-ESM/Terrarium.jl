# # [Global snowpack with idealized seasonal forcing](@id snow_global)
#
# This example demonstrates the standalone single-layer [`SnowModel`](@ref) run on a global grid,
# accelerated by GPU (if available), driven by *idealized* analytic forcing rather than reanalysis data.
# We integrate for two full annual cycles so the seasonal accumulation and melt of the snowpack is
# visible in both hemispheres.
#
# The forcing is purely illustrative: a latitude- and season-dependent near-surface air temperature, a
# constant precipitation rate that falls as snow wherever the air is below freezing, and an idealized net
# surface heat flux that cools the pack in winter and warms (melts) it in summer. The two hemispheres are
# driven out of phase.
#
# !!! note
#     This script is intended to be run directly (it is not executed during the documentation build, as a
#     two-year global integration is expensive). Magnitudes are chosen for illustration, not realism.

using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets

import RingGrids

using Oceananigans.Units: days
using Oceananigans: JLD2Writer, TimeInterval, prettytime, FieldTimeSeries, Callback, IterationInterval
using CairoMakie, GeoMakie
import DisplayAs #hide

# ## Architecture and grid
#
# Choose the GPU if one is available, otherwise fall back to the CPU.
arch = CUDA.functional() ? GPU() : CPU()
NF = Float32
@info "Setting up global snow simulation on $arch"

# We reuse the ERA5 land–sea mask at ~1° resolution to select land points, exactly as in the
# [global soil example](@ref soil_heat_global). The snowpack is a single lumped layer, so a single
# vertical level suffices.
land_sea_frac_10km = on_architecture(arch, Terrarium.load_asset(ERA5LandInvariants(), "lsm"; NF))
ring_grid = on_architecture(arch, RingGrids.FullGaussianGrid(72))
land_sea_frac_N72 = RingGrids.interpolate(ring_grid, land_sea_frac_10km)
land_mask = land_sea_frac_N72 .> 0.5
grid = ColumnRingGrid(arch, NF, UniformSpacing(N = 1), land_mask.grid, land_mask)

# The `x`-axis of the underlying `RectilinearGrid` indexes land points in ring order; we grab the
# per-column latitudes (radians) and copy them to the device so the forcing kernels can index them.
grid_lon, grid_lat = RingGrids.get_lonlats(grid.rings)
lat_masked = grid_lat[on_architecture(CPU(), land_mask)]
lat_device = on_architecture(arch, NF.(lat_masked))

# ## Idealized seasonal forcing
#
# We drive the model through the [`InputSource`](@ref) interface: `update_inputs!` is called every time
# step and (re)computes the atmospheric and surface-flux inputs analytically. All updates are written with
# `set!` so they run on the device; the closures capture only the device latitude vector and scalars.
import Terrarium: InputSource, update_inputs!, variables

struct IdealizedSeasonalForcing{NF, L} <: InputSource{NF, :idealized_seasonal_forcing}
    "Per-column latitude (radians), resident on the simulation device"
    latitude::L
    "Length of one year in seconds"
    year::NF
    "Precipitation rate applied as snowfall where the air is below freezing (m/s SWE)"
    precip_rate::NF
    "Amplitude of the idealized seasonal net surface heat flux (W/m²)"
    flux_amplitude::NF
end

# This source updates existing input variables (declared by the model) rather than adding new ones.
variables(::IdealizedSeasonalForcing) = ()

# Idealized near-surface air temperature (°C): warm at the equator, cold at the poles, with a seasonal
# swing that grows with latitude and flips sign between hemispheres.
@inline function idealized_air_temperature(φ, t, year, ::Type{NF}) where {NF}
    ω = 2NF(π) / year
    T_mean = NF(25) - NF(45) * abs(sin(φ))          # ~25 °C at equator, ~-20 °C at poles
    A_seasonal = NF(20) * abs(sin(φ))               # larger seasonal amplitude toward the poles
    return T_mean - A_seasonal * cos(ω * t) * sign(φ)  # NH coldest at t = 0, SH out of phase
end

# We update the input fields with broadcasts over the per-column latitude vector `lat`, which lives on the
# simulation device. Writing the results straight into `interior(field)` keeps everything on the device
# (using `set!` with a coordinate function would evaluate the closure over the grid nodes on the host and
# trigger scalar indexing of the device `lat` array on the GPU).
function update_inputs!(fields, src::IdealizedSeasonalForcing{NF}, clock::Terrarium.Clock, scope::Terrarium.VarPath = ()) where {NF}
    t = convert(NF, clock.time)
    lat = src.latitude
    year = src.year
    precip = src.precip_rate
    Q0 = src.flux_amplitude
    ## near-surface air temperature (reused for the snowfall test)
    T_air = (φ -> idealized_air_temperature(φ, t, year, NF)).(lat)
    interior(fields.air_temperature)[:, 1, 1] .= T_air
    ## precipitation falls as snow wherever the air is below freezing
    interior(fields.snowfall)[:, 1, 1] .= ifelse.(T_air .< zero(NF), precip, zero(NF))
    ## idealized net surface heat flux (positive upward): cooling in winter, warming (melt) in summer
    interior(fields.surface_heat_flux)[:, 1, 1] .= Q0 .* cos(2NF(π) / year * t) .* sign.(lat)
    return nothing
end

year_seconds = NF(365 * 24 * 3600)
forcing = IdealizedSeasonalForcing(lat_device, year_seconds, NF(2.0e-3 / (24 * 3600)), NF(30))

# ## Model and initial state
#
# We start from a thin, cold snowpack everywhere so that both accumulation and melt are exercised.
model = SnowModel(grid)
initializers = (snow_water_equivalent = 0.05, snow_temperature = -5.0)
integrator = initialize(model; inputs = forcing, initializers)

# Snapshot the initial snow water equivalent.
W0 = RingGrids.Field(arch, interior(integrator.state.snow_water_equivalent), grid)
fig = heatmap(W0[:, 1, 1], title = "Initial snow water equivalent (m)", colorrange = (0, 1))
DisplayAs.PNG(fig) #hide

# ## Run for two annual cycles, saving snapshots for an animation
#
# We wrap the integrator in an Oceananigans `Simulation`, attach a `JLD2Writer` that saves the snow water
# equivalent every five days, and integrate for two years with an hourly step.
output_dir = mkpath("outputs")
output_file = joinpath(output_dir, "snow_global_swe.jld2")
sim = Simulation(integrator; Δt = 3600.0, stop_time = 2 * 365days)
Terrarium.initialize!(integrator)
sim.output_writers[:swe] = JLD2Writer(
    integrator,
    (snow_water_equivalent = integrator.state.snow_water_equivalent,);
    filename = output_file,
    overwrite_existing = true,
    schedule = TimeInterval(1days),
)

# The standalone snowpack has no substrate to absorb meltwater, so in fully-melted columns the Darcy
# drainage would drive the (prognostic) snow water equivalent slightly negative; clamp it to zero after
# each step for a physical demo (the coupled `LandModel` routes this meltwater into the soil instead).
clamp_swe!(s) = (w = interior(s.model.state.snow_water_equivalent); w .= max.(w, zero(eltype(w))); nothing)
sim.callbacks[:clamp_swe] = Callback(clamp_swe!, IterationInterval(1))

# Now we can run the simulation:
run!(sim)

# We can inspect the state of `snow_water_equivalent` at the end of the simulation and see that snow has accumulated
# in the polar regions with a snow-free tropical band:
W_end = RingGrids.Field(arch, interior(integrator.state.snow_water_equivalent), grid)
fig = heatmap(W_end[:, 1, 1], title = "Snow water equivalent after two years (m)", colorrange = (0, 1))
DisplayAs.PNG(fig) #hide

# We next load the time-varying outputs using `FieldTimeSeries` and use them to create an animation:
swe_ts = FieldTimeSeries(output_file, "snow_water_equivalent")
ring0 = RingGrids.Field(swe_ts[1], grid)[:, 1]
lond = RingGrids.get_lond(ring0)
latd = RingGrids.get_latd(ring0)
let fig = Figure(size = (1200, 660))
    n_t = Observable(1)
    title = @lift string("Snow water equivalent (m) at t = ", prettytime(swe_ts.times[$n_t]))
    ax = Axis(
        fig[1, 1];
        aspect = 2,
        title,
        xticks = 0:60:360,
        yticks = -60:30:60,
        xtickformat = values -> ["$(round(Int, v))˚E" for v in values],
        ytickformat = values -> ["$(round(Int, v))˚N" for v in values],
    )
    data = @lift Matrix(RingGrids.Field(swe_ts[$n_t], grid)[:, 1])
    hm = heatmap!(ax, lond, latd, data; colorrange = (0, 1))
    Colorbar(fig[:, end + 1], hm; label = "SWE (m)")
    Makie.record(fig, joinpath("plots", "snow_global_swe.mp4"), 1:length(swe_ts.times); framerate = 10) do i
        n_t[] = i
    end
end
# ![Snow water equivalent animation](outputs/snow_global_swe.mp4)
