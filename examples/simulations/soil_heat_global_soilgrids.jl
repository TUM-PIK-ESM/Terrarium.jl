# # [Global soil heat conduction with heterogeneous soil from SoilGrids 2.0](@id soil_heat_global_soilgrids)
# This example extends the [global soil heat conduction example](@ref soil_heat_global) with a
# heterogeneous, multi-layer soil stratigraphy whose texture (sand/silt/clay fractions) is
# prescribed from the [SoilGrids 2.0](https://soilgrids.org) dataset, accessed via
# [NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl).
#
# The key new concept demonstrated here are *namespaced input variables*: each soil horizon of
# the [`SoilStratigraphy`](@ref) declares its own `sand_fraction`, `silt_fraction`, `clay_fraction`,
# and `thickness` input variables inside a namespace named after the horizon. Input sources are matched
# to these variables via namespaced names such as `:horizon1 => :sand_fraction`.

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
using NumericalEarth.DataWrangling: MetadataSet
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
# SoilGrids 2.0 provides global predictions of soil properties on a ~10 km grid for the six
# standard [GlobalSoilMap](https://www.isric.org/projects/globalsoilmap) depth intervals
# (0-5, 5-15, 15-30, 30-60, 60-100, and 100-200 cm). Sand, silt, and clay are mass fractions
# of the fine earth (< 2 mm) component reported in g/kg following the USDA particle-size
# classification (clay < 2 µm, silt 2-50 µm, sand 50-2000 µm); see the
# [GlobalSoilMap specifications (2015)](https://www.isric.org/sites/default/files/GlobalSoilMap_specifications_december_2015_2.pdf)
# for the full reference. NumericalEarth.jl downloads the data on first use (~several hundred
# MB, cached afterwards), converts the fractions to the unit interval, and returns Oceananigans
# `Field`s with the vertical axis ordered from the deepest layer (k=1, 100-200 cm) to the
# surface layer (k=6, 0-5 cm).
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
            ring_field = RingGrids.on_architecture(arch, RingGrids.FullClenshawField(interior(var_field)[:, (end - 1):-1:2, idx], input_as = Matrix))
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

soilgrids_inputs = InputSources(SoilGrids2(), grid)

fig = heatmap(RingGrids.Field(on_architecture(CPU(), soilgrids_inputs.sources[2].field), grid)[:, 1])
DisplayAs.PNG(fig) #hide

# ## Heterogeneous soil stratigraphy
# We now construct a six-layer `SoilStratigraphy` using the `SoilGridsStratigraphy` convenience constructor.
porosity = SoilPorositySURFEX(eltype(grid))
strat = SoilGridsStratigraphy(eltype(grid))
soil = SoilEnergyWaterCarbon(eltype(grid); strat)
model = SoilModel(grid; soil)

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
integrator = initialize(model, ForwardEuler(NF); inputs = soilgrids_inputs, boundary_conditions = bc, initializers = inits)

# We can verify that the SoilGrids texture has been correctly assigned to the first horizon:
sand1 = RingGrids.Field(arch, interior(integrator.state.namespaces.horizon1.sand_fraction), grid)
fig = heatmap(sand1[:, 1, 1], title = "Sand fraction of the organic horizon")
DisplayAs.PNG(fig) #hide

# Now run the simulation for 12 hours:
timestep!(integrator)
@time run!(integrator, period = Hour(12), Δt = 600.0)

# ...and plot the resulting surface temperature:
T_surface = RingGrids.Field(arch, interior(integrator.state.ground_temperature), grid)
fig = heatmap(T_surface[:, 1, 1], title = "Temperature of uppermost soil layer", colorrange = (-20, 20))
DisplayAs.PNG(fig) #hide
