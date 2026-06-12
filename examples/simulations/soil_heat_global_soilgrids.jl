# # [Global soil heat conduction with heterogeneous soil from SoilGrids 2.0](@id soil_heat_global_soilgrids)
# This example extends the [global soil heat conduction example](@ref soil_heat_global) with a
# heterogeneous, multi-layer soil stratigraphy whose texture (sand/silt/clay fractions) is
# prescribed from the [SoilGrids 2.0](https://soilgrids.org) dataset, accessed via
# [NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl).
#
# The key new concept demonstrated here are *namespaced input variables*: each soil horizon of
# the [`SoilStratigraphy`](@ref) declares its own `sand`, `silt`, `clay`, and `thickness` input
# variables inside a namespace named after the horizon. Input sources are matched to these
# variables via namespaced names such as `:organic => :sand`.

using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie

import RingGrids
import DisplayAs #hide

# NumericalEarth.jl provides the `SoilGrids2` dataset together with download and
# preprocessing machinery.
using NumericalEarth.DataWrangling: Metadatum
using NumericalEarth.SoilGrids: SoilGrids2

input_dir = "inputs" #hide
@info "Current working directory: $(pwd())" #hide

# First we check if a GPU is available and choose the architecture correspondingly.
arch = CUDA.functional() ? GPU() : CPU()
@info "Setting up simulation on $arch" #hide

# As in the [global example](@ref soil_heat_global), we load a land-sea mask at ~1° resolution
# and set up a masked [`ColumnRingGrid`](@ref) with 30 exponentially spaced soil layers.
NF = Float32
land_sea_frac = convert.(NF, dropdims(Raster(joinpath(input_dir, "era5-land_land_sea_mask_N72.nc")), dims = Ti))
land_sea_frac_field = RingGrids.FullGaussianGrid(Matrix(land_sea_frac), input_as = Matrix)
land_mask = land_sea_frac_field .> 0.5
grid = ColumnRingGrid(arch, NF, ExponentialSpacing(N = 30), land_mask.grid, land_mask)
grid_lon, grid_lat = RingGrids.get_lonlats(grid.rings) # in radians

# ## Loading soil texture from SoilGrids 2.0
# SoilGrids 2.0 provides global predictions of soil properties on a ~10 km grid for six
# depth intervals (0-5, 5-15, 15-30, 30-60, 60-100, and 100-200 cm). NumericalEarth.jl
# downloads the data on first use (~several hundred MB, cached afterwards) and returns
# Oceananigans `Field`s with the vertical axis ordered from the deepest layer (k=1,
# 100-200 cm) to the surface layer (k=6, 0-5 cm).
soilgrids = SoilGrids2()
sand3d = Field(Metadatum(:sand_fraction, dataset = soilgrids))
silt3d = Field(Metadatum(:silt_fraction, dataset = soilgrids))
clay3d = Field(Metadatum(:clay_fraction, dataset = soilgrids))

# The SoilGrids data is provided on a full Clenshaw-Curtis grid, so we can wrap it directly
# as a `RingGrids.Field` and interpolate it onto the rings of our (coarser) model grid.
# Note that NumericalEarth normalizes the latitude axis to run from south to north, whereas
# RingGrids orders rings from north to south, so we flip the latitude dimension.
function to_model_grid(field3d, k, target_grid)
    data = Array(interior(field3d, :, :, k))
    ring_field = RingGrids.FullClenshawField(reverse(data, dims = 2), input_as = Matrix)
    return RingGrids.interpolate(target_grid, ring_field)
end

# We use the topmost layer (0-5 cm) for the organic horizon and the 60-100 cm layer for
# the mineral surface horizon below it.
k_top, k_sub = 6, 2
sand_org, silt_org, clay_org = (to_model_grid(f, k_top, grid.rings) for f in (sand3d, silt3d, clay3d))
sand_sub, silt_sub, clay_sub = (to_model_grid(f, k_sub, grid.rings) for f in (sand3d, silt3d, clay3d))

# TODO: maybe move this into the package?
# The sand, silt, and clay fractions in SoilGrids are predicted independently and thus do
# not sum exactly to unity, which the [`SoilTexture`](@ref) constructor requires. We therefore
# normalize the fractions in each grid cell and fill cells without data (e.g. ocean cells
# that fall inside the land mask) with a loam-like default texture.
function normalize_texture!(sand, silt, clay; defaults = (0.4, 0.4, 0.2))
    total = sand .+ silt .+ clay
    for i in eachindex(total)
        if isnan(total[i]) || total[i] <= 0
            sand[i], silt[i], clay[i] = defaults
        else
            sand[i] /= total[i]
            silt[i] /= total[i]
            clay[i] /= total[i]
        end
    end
    return sand, silt, clay
end
normalize_texture!(sand_org, silt_org, clay_org)
normalize_texture!(sand_sub, silt_sub, clay_sub)

fig = heatmap(sand_org, title = "SoilGrids 2.0 sand fraction (0-5 cm) on the model grid")
DisplayAs.PNG(fig) #hide

# ## Heterogeneous soil stratigraphy
# We now construct a three-horizon organic/surface/bedrock (OAR) stratigraphy where the upper
# two horizons are [`PrescribedSoilHorizon`](@ref)s, i.e. their texture and thickness are
# provided as (namespaced) input variables.
strat = OARSoilStratigraphy(NF)
soil = Terrarium.SoilEnergyWaterCarbon(NF; strat)
model = SoilModel(grid; soil)

# The horizon thicknesses are also input variables; here we simply prescribe spatially
# constant values of 0.1 m for the organic horizon and 2 m for the mineral surface horizon.
# The bedrock horizon extends through the rest of the column.
organic_thickness = Field(grid, XY())
surface_thickness = Field(grid, XY())
set!(organic_thickness, NF(0.1))
set!(surface_thickness, NF(2.0))

# Each input source is assigned to its soil horizon via the namespaced `name`, e.g.
# `:organic => :sand` targets the `sand` input variable in the `organic` namespace.
inputs = InputSources(
    InputSource(grid, sand_org; name = :organic => :sand),
    InputSource(grid, silt_org; name = :organic => :silt),
    InputSource(grid, clay_org; name = :organic => :clay),
    InputSource(grid, organic_thickness; name = :organic => :thickness),
    InputSource(grid, sand_sub; name = :surface => :sand),
    InputSource(grid, silt_sub; name = :surface => :silt),
    InputSource(grid, clay_sub; name = :surface => :clay),
    InputSource(grid, surface_thickness; name = :surface => :thickness),
)

# We reuse the simple latitude-dependent climatology from the [global example](@ref soil_heat_global)
# for the initial and boundary conditions.
mean_annual_temperature(lat) = 20 - abs(40 * sin(lat))

lon_masked = grid_lon[land_mask]
lat_masked = grid_lat[land_mask]

function initial_soil_temperature(x, z)
    latᵢ = lat_masked[round(Int, x)]
    T₀ = mean_annual_temperature(latᵢ)
    T = T₀ - 0.05 * z
    return T
end

function get_temperature_bc(lon::AbstractVector, lat::AbstractVector, amplitude = 10.0)
    lon_device = on_architecture(arch, NF.(lon))
    lat_device = on_architecture(arch, NF.(lat))
    function periodic_bc(x::NF, t::NF) where {NF}
        lonₓ = lon_device[round(Int, x)]
        latₓ = lat_device[round(Int, x)]
        T₀ = mean_annual_temperature(latₓ)
        seconds_per_day = NF(24 * 3600)
        T = T₀ + NF(amplitude) * sin(2π * t / seconds_per_day - lonₓ)
        return T
    end
    return periodic_bc
end

bc = PrescribedSurfaceTemperature(:T_ub, get_temperature_bc(lon_masked, lat_masked))
inits = (temperature = initial_soil_temperature,)

# Initialize the model; the input sources are matched to the namespaced input variables of
# the prescribed soil horizons.
integrator = initialize(model, ForwardEuler(NF); inputs, boundary_conditions = bc, initializers = inits)

# We can verify that the SoilGrids texture has been correctly assigned to the organic horizon:
organic_sand = RingGrids.Field(arch, interior(integrator.state.namespaces.organic.sand), grid)
fig = heatmap(organic_sand[:, 1, 1], title = "Sand fraction of the organic horizon")
DisplayAs.PNG(fig) #hide

# Now run the simulation for 12 hours:
timestep!(integrator)
@time run!(integrator, period = Hour(12), Δt = 600.0)

# ...and plot the resulting surface temperature:
T_surface = RingGrids.Field(arch, interior(integrator.state.ground_temperature), grid)
fig = heatmap(T_surface[:, 1, 1], title = "Temperature of uppermost soil layer", colorrange = (-20, 20))
DisplayAs.PNG(fig) #hide
