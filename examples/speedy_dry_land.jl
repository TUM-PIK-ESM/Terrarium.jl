using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather as Speedy

# Temporary workaround for SpeedyWeather.jl#919
Speedy.SpeedyTransforms.FFTW.set_num_threads(1)

"""
Naive implementation of a SpeedyWeather "dry" land model that couples to a Terrarium model.
"""
struct TerrariumDryLand{NF, TMI<:Terrarium.ModelIntegrator{NF}} <: Speedy.AbstractDryLand
    "Speedy spectral grid"
    spectral_grid::Speedy.SpectralGrid

    "Initialized Terrarium model integrator"
    integrator::TMI

    function TerrariumDryLand(itnegrator::Terrarium.ModelIntegrator{NF, Arch, Grid}; spectral_grid_kwargs...) where {NF, Arch, Grid<:ColumnRingGrid}
        spectral_grid = Speedy.SpectralGrid(integrator.model.grid.rings; NF, nlayers_soil=1, spectral_grid_kwargs...)
        return new{eltype(integrator), typeof(integrator)}(spectral_grid, integrator)
    end
end

function Speedy.initialize!(
    progn_land::Speedy.PrognosticVariablesLand,  # for dispatch
    progn::Speedy.PrognosticVariables{NF},
    diagn::Speedy.DiagnosticVariables{NF},
    land::TerrariumDryLand,
    model::Speedy.PrimitiveEquation,
) where {NF}
    Tsoil = interior(land.integrator.state.temperature)[:, 1, end]
    progn_land.soil_temperature .= Tsoil .+ NF(273.15)
    Terrarium.initialize!(land.integrator)
    return nothing
end

function Speedy.timestep!(
    progn::Speedy.PrognosticVariables{NF},
    diagn::Speedy.DiagnosticVariables{NF},
    land::TerrariumDryLand,
    model::Speedy.PrimitiveEquation,
) where {NF}
    # get speedy state variables
    Tair = @view diagn.grid.temp_grid[:, end]
    # terrarium state
    state = land.integrator.state
    set!(state.inputs.air_temperature, Tair) # first set directly to avoid allocating new fields
    set!(state.inputs.air_temperature, state.inputs.air_temperature - 273.15) # then convert to celsius
    # run land forward over speedy timestep interval;
    # we use a smaller actual timestep to ensure stability
    Terrarium.run!(land.integrator, period=progn.clock.Δt, Δt=300.0)
    # Get soil temperatures
    Tsoil = state.temperature
    # Get surface temperature (last z-layer in Oceananigans grids)
    Nx, Nz = Tsoil.grid.Nx, Tsoil.grid.Nz
    Tsurf = @view Tsoil[1:Nx, 1, Nz:Nz]
    # Update speedy soil/skin temperature
    progn.land.soil_temperature .= Tsurf .+ NF(273.15)
    return nothing
end

ring_grid = RingGrids.FullGaussianGrid(24)
Nz = 30
Δz_min = 0.05 # currently the coupling is only stable with a large surface layer
grid = ColumnRingGrid(CPU(), Float32, ExponentialSpacing(; N=Nz, Δz_min), ring_grid)
# grid = ColumnGrid(CPU(), Float32, ExponentialSpacing(N=30))
# Initial conditions
soil_initializer = FieldInitializers(
    # steady-ish state initial condition for soil temperature
    temperature = (x,z) -> 0 - 0.02f0*z,
    # fully saturated soil
    saturation_water_ice = 1.0f0,
)

# Soil model with prescribed surface temperautre BC
model = SoilModel(grid, initializer=soil_initializer)
Tair_input = InputSource(eltype(grid), :air_temperature)
bcs = PrescribedSurfaceTemperature(:air_temperature)
integrator = initialize(model, ForwardEuler(eltype(grid)), boundary_conditions=bcs, inputs=InputSources(Tair_input))

# Initialize Terrarium-Speedy land model
land = TerrariumDryLand(integrator)
# Set up coupled model
land_sea_mask = Speedy.RockyPlanetMask(land.spectral_grid)
output = Speedy.NetCDFOutput(land.spectral_grid, Speedy.PrimitiveDryModel, path="outputs/")
time_stepping = Speedy.Leapfrog(land.spectral_grid, Δt_at_T31=Minute(15))
primitive_dry_coupled = Speedy.PrimitiveDryModel(land.spectral_grid; land, land_sea_mask, time_stepping, output)
# add soil temperature as output variable for Speedy simulation
Speedy.add!(primitive_dry_coupled.output, Speedy.SoilTemperatureOutput())
# initialize coupled simulation
sim_coupled = Speedy.initialize!(primitive_dry_coupled)
# Speedy.timestep!(sim)
# run it
Speedy.run!(sim_coupled, period=Month(1), output=true)

# Soil temperature in the 5th layer (~0.54 m)
Tsoil_fig = heatmap(RingGrids.Field(interior(integrator.state.temperature)[:,1,end-4], grid), title="", size=(800,400))
# Atmosphere variables
Tair_fig = heatmap(sim_coupled.diagnostic_variables.grid.temp_grid[:,8] .- 273.15, title="Air temperature", size=(800,400))
pres_fig = heatmap(exp.(sim_coupled.diagnostic_variables.grid.pres_grid), title="Surface pressure", size=(800,400))
srad_fig = heatmap(exp.(sim_coupled.diagnostic_variables.physics.surface_shortwave_down), title="Surface shortwave down", size=(800,400))
# Tskin_fig = heatmap(RingGrids.Field(interior(integrator.state.skin_temperature)[:,1,end], grid), title="", size=(800,400))
# save("plots/speedy_primitive_dry_coupled_tair_R48.png", Tair_fig, px_per_unit=1)
# save("plots/speedy_primitive_dry_coupled_tsoil_R48.png", Tsoil_fig, px_per_unit=1)

# animate surface air and soil temperatures
Speedy.animate(sim, variable="temp", coastlines=false, level=spectral_grid.nlayers, output_file="plots/speedy_terrarium_dry_air_temperature.mp4")
Speedy.animate(sim, variable="st", coastlines=false, level=1, output_file="plots/speedy_terrarium_dry_soil_temperature.mp4")

# pick a point somewhere in the mid-lattitudes
T = interior(integrator.state.temperature)[2000,1,:]
f = interior(integrator.state.liquid_water_fraction)[2000,1,:]
zs = znodes(integrator.state.temperature)
# Plot temperature and liquid fraction profiles in upper 15 layers
Makie.scatterlines(T[end-15:end], zs[end-15:end])
Makie.scatterlines(f, zs)