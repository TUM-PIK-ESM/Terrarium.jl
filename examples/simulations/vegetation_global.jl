# # [Simulating global vegetation with prescribed LAI](@id vegetation_prescribed_global)
# This example sets up a global [`VegetationModel`](@ref) driven by the [`PrescribedPhenology`](@ref)
# scheme, using the ERA5-Land leaf area index (LAI) climatology for high vegetation (`lai_hv`) as a
# annually repeating, time-varying input. The phenology factor `ϕ` is derived from the prescribed LAI
# and used by the vegetation carbon-cycle processes.

using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Rasters: X, Y, Ti

using CairoMakie, GeoMakie

import RingGrids

import DisplayAs #hide

# Run on GPU if available, otherwise CPU.
arch = CUDA.functional() ? GPU() : CPU()
NF = Float32
@info "Setting up simulation on $arch" #hide

# ## Grid and land mask
# As in the [global soil heat example](@ref soil_heat_global), we build a masked
# [`ColumnRingGrid`](@ref) at ~1° resolution (72 Gaussian rings), keeping only grid points with
# >50% land cover. The land-sea mask is loaded from the ERA5-Land invariants asset and regridded
# onto the model grid.
land_sea_frac_native = on_architecture(arch, Terrarium.load_asset(ERA5LandInvariants(), "lsm"; NF))
model_rings = on_architecture(arch, RingGrids.FullGaussianGrid(72))
land_sea_frac = RingGrids.interpolate(model_rings, land_sea_frac_native)
land_mask = land_sea_frac .> 0.5
grid = ColumnRingGrid(arch, NF, ExponentialSpacing(N = 10), land_mask.grid, land_mask)

# ## Prescribed LAI climatology
# The `ERA5LandLeafAreaIndex` asset holds a daily LAI climatology (1980–2010) at the native 0.1°
# resolution. Note that, due to the high resolution, it may take a few minutes to download the data.
lai_asset = ERA5LandLeafAreaIndex()
lai_asset_path = Terrarium.get_asset(lai_asset)

# We next load the asset into a `RasterStack` with `lazy = true` to avoid loading the full dataset.
# cycle = true` makes the LAI climatology repeat every year over the whole simulation.
lai_raster = RasterStack(lai_asset_path, lazy = true)
lai_highveg = replace_missing(lai_raster[:lai_hv], NaN32)
lai_input = InputSource(grid, lai_highveg; name = :leaf_area_index, cycle = true)

# ## Vegetation model
# We build a [`VegetationCarbonCycle`](@ref) with the [`PrescribedPhenology`](@ref) scheme. Prescribed
# phenology imposes LAI externally rather than deriving it from a prognostic carbon pool, so we disable
# the prognostic vegetation-fraction dynamics (`vegetation_dynamics = nothing`). The carbon pool
# `carbon_vegetation` is still needed by the respiration parameterization, so we initialize it to a
# constant; note that with prescribed LAI this pool is not consistently coupled to the imposed leaf area.
vegetation = VegetationCarbonCycle(
    NF;
    phenology = PrescribedPhenology(NF; max_leaf_area_index = 6),
    vegetation_dynamics = nothing,
)
model = VegetationModel(grid; vegetation)

# Initialize the model with the LAI input source and a uniform initial vegetation carbon pool.
integrator = initialize(model; inputs = lai_input, initializers = (carbon_vegetation = NF(0.5),))

# Helper to plot a horizontal field on the ring grid.
plot_field(field; kwargs...) = heatmap(RingGrids.Field(arch, interior(field), grid)[:, 1, 1]; kwargs...)

# The initial (mid-January) prescribed LAI and derived phenology factor:
fig = plot_field(integrator.state.leaf_area_index, title = "Leaf area index (Jan)", colorrange = (0, 6))
DisplayAs.PNG(fig) #hide

# ## Run through the growing season
# We advance the simulation for six months with a daily timestep. The cyclic LAI input drives the
# seasonal evolution of the phenology factor.
@time run!(integrator, period = Day(180), Δt = 24 * 3600.0)

# The mid-summer LAI after half a year:
fig = plot_field(integrator.state.leaf_area_index, title = "Leaf area index (Jul)", colorrange = (0, 6))
DisplayAs.PNG(fig) #hide
