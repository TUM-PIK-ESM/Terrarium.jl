# # Global sensitivity analysis with Reactant + Enzyme
#
# This example computes a *global sensitivity map*: the sensitivity of the global mean surface
# soil temperature after one simulated day with respect to the initial soil state at **every**
# land grid point and soil layer. Reverse-mode AD delivers this entire map in a *single* gradient
# evaluation — the cost of the reverse pass is independent of the number of inputs, which is what
# makes global (all-points-at-once) sensitivity analysis feasible.
#
# We follow the same Reactant + Enzyme recipe as the
# [Reactant differentiation example](@ref): the model is built on `ReactantState()`, the stepping
# loop [`run_timesteps!`](@ref) is traced with a checkpointing scheme, and `Enzyme.autodiff`
# differentiates through the whole compiled rollout. The physical setup mirrors the ERA5-forced
# global soil heat example (`examples/simulations/soil_heat_global_era5.jl`).

using Terrarium
using Reactant, CUDA   # CUDA is required by Reactant's kernel integration, even on CPU
using Enzyme

using Rasters, NCDatasets
using CairoMakie
using Statistics: quantile

import RingGrids

# ## Grid and forcing
#
# We build the same ~1° global land grid as the ERA5 example: all grid points with more than 50%
# land become independent soil columns. Reactant currently requires uniform vertical spacing (see
# the [Reactant page](@ref)), so we use `UniformSpacing` instead of `ExponentialSpacing`.

## Load land-sea mask at ~1° resolution
land_sea_frac = convert.(Float32, dropdims(Raster("inputs/era5-land_land_sea_mask_N72.nc"), dims = Ti))
land_sea_frac_field = RingGrids.FullGaussianField(Matrix(land_sea_frac), input_as = Matrix)
land_mask = land_sea_frac_field .> 0.5 # select only grid points with > 50% land

Nz = 30 # number of soil layers
grid = ColumnRingGrid(ReactantState(), Float32, UniformSpacing(Δz = 0.2f0, N = Nz), land_mask)

# Load ERA-5 2 meter air temperature at ~1° resolution as the surface forcing. The compiled
# stepping loop cannot yet update time-dependent `FieldTimeSeries` inputs from the host, so we
# take the first time slice as a constant-in-time forcing field (a `FieldInputSource`, as in the
# hybrid modeling example).

# The full-year hourly file is large (~1.5 GB, 8760 time steps); load it lazily and materialize
# only the first time slice, which we use as a constant-in-time forcing.
Tair_raster = Raster("inputs/era5_land_2m_temperature_2023_N72.nc"; lazy = true)
Tair0 = convert.(Float32, replace_missing(Tair_raster[Ti(1)], NaN)) .- 273.15f0
Tair_field = RingGrids.FullGaussianField(Matrix(Tair0), input_as = Matrix)
Tair_forcing = InputSource(grid, Tair_field, name = :Tair, units = u"°C")

# Masked surface temperatures ordered by the grid's active (land) columns, for the initializer
# below (column index `x` runs over the land columns).
Tsurf_0 = Tair_field[findall(land_mask)]

# ## Model setup
#
# As everywhere with Reactant in Terrarium, the model itself is unchanged — only the grid's
# architecture differs. Initialization runs eagerly and is transferred to the device.

model = SoilModel(grid, timestepper = ForwardEuler(eltype(grid)))
boundary_conditions = PrescribedSurfaceTemperature(:Tair)
initializers = (
    ## steady-ish state initial condition for temperature
    temperature = (x, z) -> Tsurf_0[Int(round(x)) + 1] - 0.02f0 * z,
    saturation_water_ice = 1.0f0,
)
inputs = InputSources(Tair_forcing)
integrator = initialize(model; inputs, initializers, boundary_conditions)

# ## The objective and its gradient
#
# The scalar objective is the global mean surface (top-layer) soil temperature after `nsteps`.
# [`run_timesteps!`](@ref) is the stepping loop underlying `run!`; on `ReactantState` it traces
# the loop and takes a `checkpointing` scheme that bounds the memory of the reverse pass.

function mean_surface_temperature(integrator, Δt, nsteps, checkpointing)
    run_timesteps!(integrator, Δt, nsteps, checkpointing)
    Tsurf = interior(integrator.state.temperature)[:, :, end]
    return sum(Tsurf) / length(Tsurf)
end

# `Enzyme.autodiff` in reverse mode computes the objective together with its gradient. The
# integrator is passed as `Duplicated`, so its shadow `dintegrator` accumulates the sensitivity
# of the objective with respect to every state variable — in particular the prognostic internal
# energy `U₀` of every column and layer.

function grad_mean_surface_temperature!(integrator, dintegrator, Δt, nsteps, checkpointing)
    _, value = Enzyme.autodiff(
        Enzyme.set_strong_zero(Enzyme.ReverseWithPrimal),
        mean_surface_temperature, Enzyme.Active,
        Enzyme.Duplicated(integrator, dintegrator),
        Enzyme.Const(Δt),
        Enzyme.Const(nsteps),
        Enzyme.Const(checkpointing),
    )
    return value
end

# We allocate the shadow with `make_zero` and compile the whole gradient with Reactant. As with a
# forward `run!`, the first call compiles (this can take a few minutes); afterwards it is fast.

dintegrator = Enzyme.make_zero(integrator)
Δt = 600.0f0
nsteps = 144                                       # one simulated day
checkpointing = Reactant.Periodic(isqrt(nsteps))   # ≈ √n checkpoints

compiled_grad! = @compile raise = true raise_first = true sync = true grad_mean_surface_temperature!(
    integrator, dintegrator, Δt, nsteps, checkpointing
)
T̄_surf = compiled_grad!(integrator, dintegrator, Δt, nsteps, checkpointing)
println("Global mean surface soil temperature after 1 day: $(Reactant.to_number(T̄_surf)) °C")

# ## The global sensitivity map
#
# A single reverse pass gave us ``\partial \bar{T}_\text{surf} / \partial U_0`` for all land
# columns and soil layers at once. We move the shadow to the CPU and scatter it back onto the
# ring grid for plotting.

dU = Array(interior(dintegrator.state.internal_energy))   # (Ncolumns, 1, Nz)
cpu_grid = on_architecture(CPU(), grid)
## scatter all layers back onto the rings (columns × layers), then take the top (surface) layer as
## a horizontal field for plotting.
dU_ring = RingGrids.Field(dU[:, 1, :], cpu_grid)
dU_surface = dU_ring[:, end]                              # top-layer sensitivity, a horizontal field

# The sensitivity is positive almost everywhere (more initial internal energy ⇒ warmer soil ⇒
# higher global mean), but a small number of columns sitting near the freezing point have an
# outsized response: there latent heat of fusion makes the effective heat capacity very large and
# ``\partial T / \partial U`` spikes. Those outliers are physical, but on a linear color scale
# they would pin every other column to one end of the colorbar. We therefore clip the color range
# to the 2nd–98th percentile so the spatial structure of the bulk is visible.
lo, hi = quantile(filter(isfinite, Array(dU_surface)), (0.02, 0.98))

mkpath("plots")
fig = heatmap(
    dU_surface;
    title = "Sensitivity of global mean surface soil temperature to initial internal energy",
    size = (900, 450),
    colorrange = (lo, hi),
    colormap = :viridis,
    highclip = :magenta,
)
Makie.save("plots/global_sensitivity_initialU_era5land.png", fig)
fig

# The same gradient also tells us how deep into the soil a single day of surface forcing
# "reaches": we average the magnitude of the sensitivity over all land columns by depth. The
# sensitivity falls off by many orders of magnitude within the first metre — a single day of
# surface forcing barely perturbs the deep soil — so we plot the magnitude on a log axis.

zs = znodes(integrator.state.temperature)
dU_depth = vec(sum(abs, dU[:, 1, :], dims = 1)) ./ size(dU, 1)
## floor at a tiny fraction of the maximum so the log axis is well-defined without flattening the
## many-orders-of-magnitude decay with depth (values span ~1e-12 at the surface to ~1e-26 deep).
floor_value = maximum(dU_depth) * 1.0f-20
dU_depth_plot = max.(dU_depth, floor_value)

f = Makie.Figure()
Makie.Axis(
    f[1, 1], ylabel = "Soil depth / m", xlabel = "mean |∂T̄_surf/∂U₀|",
    xscale = log10,
)
Makie.scatterlines!(f[1, 1], dU_depth_plot, zs)
Makie.save("plots/global_sensitivity_depth_profile.png", f)
f

# Sensitivities with respect to *model parameters* (rather than initial conditions) are shown in
# the [differentiation example](@ref); with Reactant the same `Duplicated(integrator, dintegrator)`
# shadow carries them too, as demonstrated for the neural-network weights in the hybrid modeling
# example (`examples/hybrid_models/neural_snow_melt.jl`).
