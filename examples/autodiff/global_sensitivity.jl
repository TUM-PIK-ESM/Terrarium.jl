using Terrarium

using CUDA
using Rasters, NCDatasets

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather

# run on GPU if available
arch = CPU() #CUDA.functional() ? GPU() : CPU()

# Load land-sea mask at ~1° resolution
land_sea_frac = convert.(Float32, dropdims(Raster("inputs/era5-land_land_sea_mask_N72.nc"), dims = Ti))
land_sea_frac_field = RingGrids.FullGaussianField(Matrix(land_sea_frac), input_as = Matrix)
heatmap(land_sea_frac_field)

# Load ERA-5 2 meter air temperature at ~1° resolution
Tair_raster = Raster("inputs/external/era5-land/2m_temperature/era5_land_2m_temperature_2023_N72.nc")
Tair_raster = convert.(Float32, replace_missing(Tair_raster, NaN)) .- 273.15f0
# heatmap(Tair_raster[:,:,1])

# Set up grids
land_mask = land_sea_frac_field .> 0.5 # select only grid points with > 50% land
Nz = 30 # number of soil layers
grid = ColumnRingGrid(arch, ExponentialSpacing(N = Nz), land_mask)

# Construct input sources
Tair_forcing = InputSource(grid, rebuild(Tair_raster, name = :Tair))
Tsurf_0 = Tair_raster[Ti(1)][findall(land_mask)]

model = SoilModel(grid)
boundary_conditions = PrescribedSurfaceTemperature(:Tair)
# Initial conditions
initializers = (
    # steady-ish state initial condition for temperature
    temperature = (x, z) -> Tsurf_0[Int(round(x)) + 1] - 0.02 * z,
    # dry soil
    saturation_water_ice = 1.0,
)
integrator = initialize(model, ForwardEuler(eltype(grid)), Tair_forcing; initializers, boundary_conditions)
@time timestep!(integrator)
@time run!(integrator, period = Hour(1), Δt = 120.0)

# plot heatmap of soil temperature at the surface
heatmap(RingGrids.Field(integrator.state.temperature, grid)[:, end])
heatmap(RingGrids.Field(integrator.state.liquid_water_fraction, grid)[:, end])

using Enzyme
using Statistics

function mean_surface_temperature(clock, model, inputs, state, inits, timestepper)
    integrator = Terrarium.ModelIntegrator(clock, model, inputs, state, inits, timestepper)
    run!(integrator, steps = 1)
    return mean(interior(integrator.state.temperature)[:, :, end])
end

dmodel = Ref(make_zero(model))
dintegrator = make_zero(integrator)
dstate = dintegrator.state
dinputs = dintegrator.inputs

autodiff(
    set_runtime_activity(Enzyme.Reverse),
    mean_surface_temperature,
    Active,
    Const(integrator.clock),
    MixedDuplicated(model, dmodel),
    Duplicated(integrator.inputs, dinputs),
    Duplicated(integrator.state, dstate),
    Const(integrator.initializers),
    Duplicated(integrator.timestepper, make_zero(integrator.timestepper))
)

@show dmodel.x.soil.energy.thermal_properties.conductivities
# compute sensitivity of T to the parameter of interest as:
# ∂f∂p⁻¹⋅∂f∂T
dTdp = dmodel.x.soil.energy.thermal_properties.conductivities.ice / dstate.temperature

zs = znodes(integrator.state.temperature)
Makie.with_theme(fontsize = 18) do
    fig = heatmap(RingGrids.Field(Field(dstate.temperature), grid)[:, end], title = "Sensitivity of global mean soil surface temperature to initial condition", size = (800, 400))
    Makie.save("plots/global_sensitivity_initialT_era5land.png", fig)
    fig
end
