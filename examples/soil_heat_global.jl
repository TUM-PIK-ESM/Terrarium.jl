using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather

const num_soil_layers = 30

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
grid = ColumnRingGrid(GPU(), ExponentialSpacing(N=num_soil_layers), land_mask)

# Construct input sources
forcings = InputSource(grid, rebuild(Tair_raster, name=:Tair))
# Calculate mean annual air temperature assuming hourly time resolution
Tsurf_avg = aggregate(mean, Tair_raster, (Ti(365*24), X(1), Y(1)))
Rasters.rplot(Tsurf_avg)

# Initial conditions
Tsurf_initial = Tsurf_avg[findall(land_mask)]
thermal_constant = sqrt(π / (1e-5 * 24*3600*365)) # with arbitrary thermal diffusivity
initializer = FieldInitializers(
    # steady-ish state initial condition calculated from mean annual air temperature;
    # consider moving to built-in initializer type?
    temperature = (x,z) -> Tsurf_initial[Int(round(x))+1] - 0.025*z + 20*exp(z*thermal_constant)*sin(z*thermal_constant),
    # dry soil
    pore_water_ice_saturation = 1.0,
)
T_ub = PrescribedValue(:temperature, Input(:Tair, units=u"°C"))
boundary_conditions = SoilBoundaryConditions(grid, top=T_ub)
model = SoilModel(grid; initializer, boundary_conditions)
state = initialize(model, forcings)
# advance one timestep with Δt = 15 minutes
@time timestep!(state, Minute(15))
# run multiple timesteps over a given time period
@time run!(state, period=Day(366), Δt=Minute(15))

function zonal_nanmean(field::RingGrids.AbstractField)
    # extend to any non-horizontal dimensions of grid (e.g. vertical or time)
    ks = size(field)[2:end]

    # determine type T after division with integer (happening in mean)
    T = Base.promote_op(/, eltype(field), Int64)
    m = zeros(T, RingGrids.get_nlat(field), ks...)

    rings = RingGrids.eachring(field.grid)
    for k in RingGrids.eachlayer(field)
        for (j, ring) in enumerate(rings)
            m[j, k] = mean(skipmissing(replace(field[ring, k], NaN => missing)))
        end
    end

    return m
end

function plot_zonal_heatmap!(fig, tempobs::Observable, t_obs::Observable; max_depth=10.0)
    Tsoil_zonal = @lift reduce(hcat, [zonal_nanmean(RingGrids.Field(CPU(), view($tempobs, :, :, l), grid)) for l in 1:num_soil_layers])
    zone_mask = @lift any(isfinite.($Tsoil_zonal), dims=2)[:,1]
    zones = @lift RingGrids.get_latd(grid.rings)[$zone_mask]
    zs = Array(znodes(tempobs[]))
    zmask = zs .> -max_depth
    Tsoil = @lift $Tsoil_zonal[$zone_mask, zmask]
    ax = Axis(fig, xlabel="Latitude / °N", ylabel="Elevation / m", title=@lift("$($t_obs) days since January 1st"))
    hm = heatmap!(ax, zones, zs[zmask], Tsoil)
    return ax, hm
end

# restart simulation and record temperature state to make animation
state = initialize(model, forcings)
# Tsoil_offset = CenterField(on_architecture(CPU(), grid.grid))
# Tsurf_avg_vec = Tsurf_avg[findall(land_mask)]
Tsoil = on_architecture(CPU(), state.state.temperature)
Tsoil_initial = deepcopy(Tsoil)
ΔTsoil = Field(Tsoil - Tsoil_initial)
temperature_observable = Observable(ΔTsoil)
time_observable = Observable(state.clock.time)
fig = Figure();
ax, hm = plot_zonal_heatmap!(fig[1,1], temperature_observable, time_observable);
hm.colorrange[] = (-20.0, 20.0)
hm.colormap[] = :bwr
Colorbar(fig[1,2], hm, label="Temperature / °C");

CairoMakie.record(fig, "plots/zonal_soil_temperature.mp4", 1:364*4; framerate=10) do i
    run!(state, period=Hour(6), Δt=900.0)
    Tsoil = on_architecture(CPU(), state.state.temperature)
    set!(ΔTsoil, Tsoil - Tsoil_initial)
    temperature_observable[] = ΔTsoil
    time_observable[] = state.clock.time / (24*3600)
end

using Oceananigans.OutputWriters: JLD2Writer, AveragedTimeInterval
using Oceananigans.Utils: days

# set up and run an Oceananigans Simulation
state = initialize(model, forcings)
sim = Simulation(state; Δt=900.0, stop_iteration=100)
# sim.output_writers[:temperature] = JLD2Writer(state, get_fields(state, :temperature); filename="output/soil_model_temperature.jld2", schedule=AveragedTimeInterval(1days))
run!(sim)
