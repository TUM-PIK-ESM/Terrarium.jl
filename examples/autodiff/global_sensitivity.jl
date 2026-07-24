# # Global sensitivity analysis with Reactant + Enzyme
#
# This example computes *global sensitivity maps*: Both parameter sensitivty and state sensitivty. 
# We start by computing the sensitivity of the global mean surface
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
# As everywhere with Reactant in Terrarium, the model construction is standard — only the grid's
# architecture differs. We give the soil a **loam** texture (sand < 1): it is a representative
# composition and, in particular, it makes the non-quartz mineral thermal conductivity influence the
# bulk conductivity, which the parameter sensitivity below relies on. The same soil is used
# throughout. 

stratigraphy = HomogeneousSoilStratigraphy(
    eltype(grid), texture = SoilTexture(eltype(grid); sand = 0.4f0, clay = 0.3f0, silt = 0.3f0)
)
model = SoilModel(
    grid; soil = SoilEnergyWaterCarbon(eltype(grid); strat = stratigraphy),
    timestepper = ForwardEuler(eltype(grid))
)
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
# columns and soil layers at once. We move the shadow to the CPU and plot it. 

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

# ## Parameter sensitivity: soil mineral thermal conductivity
#
# The maps above are sensitivities to the *initial state*. We now compute the sensitivity to a
# *model parameter* — the soil mineral (non-quartz) thermal conductivity ``\kappa_\text{mineral}`` —
# both as a global number and as a per-grid-point map. A scalar physical parameter needs a little
# extra care under Reactant: on a `ReactantState` grid the state arrays and the clock are traced, but
# the model's scalar parameters are baked into the compiled program as constants, so Enzyme cannot
# see them by default. (Array-valued parameters such as neural-network weights are traced as device
# arrays and *are* differentiated directly, as in the hybrid modeling example — the asymmetry is that
# arrays decouple from the grid's number type while a scalar `::NF` field does not.)
#
# To make ``\kappa_\text{mineral}`` a differentiable input we promote the soil thermal conductivities
# to traced device scalars with `Reactant.to_rarray(...; track_numbers = true)`. 
#
# Two physical points are worth noting. First, the mineral endpoint only enters the bulk conductivity
# through the *non-quartz* mineral fraction,
# ```math
# \kappa_\text{mineral grains} = \kappa_\text{quartz}^{\,q}\,\kappa_\text{mineral}^{\,1 - q},
# ```
# where ``q`` is the sand (quartz) fraction: with a pure-sand texture (``q = 1``) the
# ``\kappa_\text{mineral}`` term would drop out entirely and the sensitivity would be *exactly* zero —
# the loam soil built above (``q = 0.4``) is what makes it nonzero. Second, we differentiate the
# **column-mean** soil temperature rather than the surface temperature: over one day the surface layer
# is held by the atmospheric boundary condition and is essentially insensitive to soil conductivity,
# whereas the conductivity controls how heat is redistributed *through* the column (the sensitivity
# peaks around 0.5 m depth and again in the deep soil).

## Build an integrator whose mineral conductivity `κ` is a *traced* device scalar (so Enzyme,
## and a compiled forward evaluation, see it as a differentiable/variable input rather than a baked
## constant). Everything else stays `Float32`; only the conductivities' number type is promoted.
function traced_integrator(κ_mineral)
    conductivities = Reactant.to_rarray(
        SoilThermalConductivities(eltype(grid); mineral = κ_mineral); track_numbers = true
    )
    thermal_properties = SoilThermalProperties(eltype(grid); conductivities)
    energy = SoilThermodynamics(eltype(grid); thermal_properties)
    soil = SoilEnergyWaterCarbon(eltype(grid); strat = stratigraphy, energy)
    model = SoilModel(grid; soil, timestepper = ForwardEuler(eltype(grid)))
    return initialize(model; inputs, initializers, boundary_conditions)
end

## objective: column-mean soil temperature (averaged over every land column and layer)
function mean_column_temperature(integrator, Δt, nsteps, checkpointing)
    run_timesteps!(integrator, Δt, nsteps, checkpointing)
    T = interior(integrator.state.temperature)
    return sum(T) / length(T)
end

function grad_mean_column_temperature!(integrator, dintegrator, Δt, nsteps, checkpointing)
    _, value = Enzyme.autodiff(
        Enzyme.set_strong_zero(Enzyme.ReverseWithPrimal),
        mean_column_temperature, Enzyme.Active,
        Enzyme.Duplicated(integrator, dintegrator),
        Enzyme.Const(Δt), Enzyme.Const(nsteps), Enzyme.Const(checkpointing),
    )
    return value
end

# ### Global sensitivity
#
# One reverse pass accumulates the parameter sensitivity into the *model* shadow (exactly as the
# neural-network weights are in the hybrid modeling example); we read it straight from the
# conductivities shadow.

κ₀ = eltype(grid)(2)                       # default mineral thermal conductivity, 2 W m⁻¹ K⁻¹
param_integrator = traced_integrator(κ₀)
param_dintegrator = Enzyme.make_zero(param_integrator)
compiled_param_grad! = @compile raise = true raise_first = true sync = true grad_mean_column_temperature!(
    param_integrator, param_dintegrator, Δt, nsteps, checkpointing
)
compiled_param_grad!(param_integrator, param_dintegrator, Δt, nsteps, checkpointing)

∂T̄_∂κ_mineral = Reactant.to_number(
    param_dintegrator.model.soil.energy.thermal_properties.conductivities.mineral
)
println("∂T̄_col/∂κ_mineral = $(∂T̄_∂κ_mineral) K / (W m⁻¹ K⁻¹)")

# ### Per-grid-point sensitivity map
#
# The global number is a spatial average; underneath it lies the sensitivity of *each* column's mean
# temperature, ``\partial \bar{T}_\text{col}(x)/\partial\kappa_\text{mineral}``. Because
# ``\kappa_\text{mineral}`` is a single scalar and the output is one value per column, *forward*-mode
# AD is the natural tool — one forward pass propagates a unit ``\kappa_\text{mineral}`` tangent to the
# whole temperature field at once (reverse mode would need one pass per column). We seed a shadow
# integrator that is zero everywhere except a `1` on the mineral conductivity, differentiate the
# rollout with `Enzyme.autodiff(Forward, …)`, and read the temperature tangent
# ``\partial T/\partial\kappa_\text{mineral}`` straight out of that shadow — the forward-mode mirror of
# the reverse pass above.

function step_rollout!(integrator, Δt, nsteps)
    run_timesteps!(integrator, Δt, nsteps)
    return nothing
end

function temperature_tangent!(integrator, seed, Δt, nsteps)
    Enzyme.autodiff(
        Enzyme.Forward, step_rollout!,
        Enzyme.Duplicated(integrator, seed),
        Enzyme.Const(Δt), Enzyme.Const(nsteps),
    )
    return nothing
end

## a fresh integrator (the reverse pass above advanced `param_integrator` past t = 0) plus a shadow
## seeded with a unit tangent on κ_mineral and zero everywhere else. The seed reuses the shadow's own
## (traced) conductivities type, so `Terrarium.setproperties` only changes the mineral *value* 0 → 1.
fwd_integrator = traced_integrator(κ₀)
seed = Enzyme.make_zero(fwd_integrator)
seed_conductivities = Terrarium.setproperties(
    seed.model.soil.energy.thermal_properties.conductivities,
    (mineral = one(seed.model.soil.energy.thermal_properties.conductivities.mineral),),
)
seed = Terrarium.setproperties(seed, (model = Terrarium.setproperties(seed.model,
    (soil = Terrarium.setproperties(seed.model.soil,
        (energy = Terrarium.setproperties(seed.model.soil.energy,
            (thermal_properties = Terrarium.setproperties(seed.model.soil.energy.thermal_properties,
                (conductivities = seed_conductivities,)),)),)),)),))

compiled_tangent! = @compile raise = true raise_first = true sync = true temperature_tangent!(
    fwd_integrator, seed, Δt, nsteps
)
compiled_tangent!(fwd_integrator, seed, Δt, nsteps)
dT_field = Array(interior(seed.state.temperature))   # ∂T/∂κ_mineral, (Ncolumns, 1, Nz)

## column-mean over layers → one sensitivity per column
Nz_soil = size(dT_field, 3)
dκ_column = vec(sum(dT_field[:, 1, :], dims = 2)) ./ Nz_soil
println("map area-average = $(sum(dκ_column) / length(dκ_column))  (matches the global scalar above)")

# Bring the per-column sensitivities back onto the ring grid and plot (broadcasting the
# column-mean across layers so the top slice `[:, end]` recovers it as a horizontal field, mirroring
# the initial-condition map).
dκ_ring = RingGrids.Field(repeat(dκ_column, 1, Nz_soil), cpu_grid)
dκ_map = dκ_ring[:, end]
lo_κ, hi_κ = quantile(filter(isfinite, Array(dκ_map)), (0.02, 0.98))
fig_κ = heatmap(
    dκ_map;
    title = "Sensitivity of column-mean soil temperature to mineral thermal conductivity",
    size = (900, 450),
    colorrange = (lo_κ, hi_κ),
    colormap = :viridis,
    highclip = :magenta,
)
Makie.save("plots/global_sensitivity_mineral_conductivity_map.png", fig_κ)
fig_κ

dκ = param_dintegrator.model.soil.energy.thermal_properties.conductivities
∂T̄_∂κ_mineral = Reactant.to_number(dκ.mineral)
println("∂T̄_surf/∂κ_mineral = $(∂T̄_∂κ_mineral) K / (W m⁻¹ K⁻¹)")
