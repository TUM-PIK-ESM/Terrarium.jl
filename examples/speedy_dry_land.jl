using Terrarium

using CUDA
using Dates
using Rasters, NCDatasets
using Statistics

using CairoMakie, GeoMakie

import RingGrids
import SpeedyWeather as Speedy

"""
Naive implementation of a SpeedyWeather "dry" land model that couples to a Terrarium model.
"""
struct TerrariumDryLand{NF, TM<:ModelState{NF}} <: Speedy.AbstractDryLand
    "Speedy spectral grid"
    spectral_grid::Speedy.SpectralGrid

    "Initialized Terrarium model"
    model::TM

    function TerrariumDryLand(model::ModelState{NF, Arch, Grid}; spectral_grid_kwargs...) where {NF, Arch, Grid<:ColumnRingGrid}
        spectral_grid = Speedy.SpectralGrid(model.grid.rings; NF, nlayers_soil=1, spectral_grid_kwargs...)
        return new{eltype(model), typeof(model)}(spectral_grid, model)
    end
end

function Speedy.initialize!(
    progn_land::Speedy.PrognosticVariablesLand,  # for dispatch
    progn::Speedy.PrognosticVariables{NF},
    diagn::Speedy.DiagnosticVariables{NF},
    land::TerrariumDryLand,
    model::Speedy.PrimitiveEquation,
) where {NF}
    Tsoil = interior(land.model.state.temperature)[:, 1, end]
    progn_land.soil_temperature .= Tsoil .+ NF(273.15)
    Terrarium.initialize!(land.model)
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
    pres = diagn.grid.pres_grid
    SwIn = diagn.physics.surface_shortwave_down
    LwIn = diagn.physics.surface_longwave_down
    SwOut = diagn.physics.surface_shortwave_up
    LwOut = diagn.physics.surface_longwave_up
    H = diagn.physics.sensible_heat_flux
    # update Terrarium input fields
    inputs = land.model.state.inputs
    set!(inputs.air_temperature, Tair) # first set directly to avoid allocating new fields
    set!(inputs.air_temperature, inputs.air_temperature - 273.15) # then convert to celsius
    # set!(inputs.air_pressure, pres) # similar with pressure
    # set!(inputs.air_pressure, exp(inputs.air_pressure)) # take exp to get pressure in Pa
    # set!(inputs.surface_shortwave_down, SwIn)
    # set!(inputs.surface_longwave_down, LwIn)
    # set!(inputs.surface_shortwave_up, SwOut)
    # set!(inputs.surface_longwave_up, LwOut)
    # set!(inputs.sensible_heat_flux, H)
    # do land timestep
    Terrarium.timestep!(land.model, Dates.second(progn.clock.Δt))
    # SEB coupling:
    # Ts = land.model.state.skin_temperature
    # progn.land.soil_temperature .= Ts .+ NF(273.15)
    # progn.land.sensible_heat_flux .= land.model.state.sensible_heat_flux
    # Soil surface temperature coupling:
    Tsoil = land.model.state.temperature
    Nx, Nz = Tsoil.grid.Nx, Tsoil.grid.Nz
    Tsurf = @view Tsoil[1:Nx, 1, Nz:Nz]
    progn.land.soil_temperature .= Tsurf .+ NF(273.15)
    return nothing
end

ring_grid = RingGrids.FullGaussianGrid(24)
Nz = 30
Δz_min = 0.1 # currently the coupling is only stable with a large surface layer
grid = ColumnRingGrid(CPU(), Float32, ExponentialSpacing(;N=Nz, Δz_min), ring_grid)
# grid = ColumnGrid(CPU(), Float32, ExponentialSpacing(N=30))
# Initial conditions
soil_initializer = FieldInitializers(
    # steady-ish state initial condition for soil temperature
    temperature = (x,z) -> 10 - 0.02f0*z,
    # fully saturated soil
    saturation_water_ice = 1.0f0,
)

# Coupled model with prescribed surface energy balance
# currently not working :(
# surface_energy_balance = SurfaceEnergyBalance(
#     eltype(grid);
#     radiative_fluxes = PrescribedRadiativeFluxes(),
#     turbulent_fluxes = PrescribedTurbulentFluxes(),
# )
# soil_model = SoilModel(grid, initializer=soil_initializer)
# model = CoupledSoilAtmosphereModel(soil_model; surface_energy_balance)
# # Quick check to make sure everything works OK
# land_state = initialize(model)
# set!(land_state.state.surface_shortwave_down, 200.0)
# set!(land_state.state.surface_longwave_down, 100.0)
# set!(land_state.state.air_temperature, 11.0)
# set!(land_state.state.air_pressure, 101_325)
# set!(land_state.state.specific_humidity, 0.001)
# timestep!(land_state)
# run!(land_state, period=Day(1))
# @assert all(isfinite.(land_state.state.skin_temperature))
# @assert all(isfinite.(land_state.state.latent_heat_flux))
# @assert all(isfinite.(land_state.state.temperature))

# Soil model with prescribed surface temperautre BC
soil_bcs = SoilBoundaryConditions(
    eltype(grid);
    top=PrescribedValue(:temperature, Input(:air_temperature, units=u"°C")),
    bottom=DefaultBoundaryConditions()
)
model = SoilModel(grid, initializer=soil_initializer, boundary_conditions=soil_bcs)

# Initialize Terrarium model
land_state = initialize(model)
land = TerrariumDryLand(land_state)

# Set up coupled model
land_sea_mask = Speedy.RockyPlanetMask(land.spectral_grid)
# surface_heat_flux = Speedy.SurfaceHeatFlux(land=Speedy.PrescribedLandHeatFlux())
output = Speedy.NetCDFOutput(land.spectral_grid, Speedy.PrimitiveDryModel, path="outputs/")
primitive_dry_coupled = Speedy.PrimitiveDryModel(land.spectral_grid; land, land_sea_mask, output)
# add soil temperature as output variable for Speedy simulation
Speedy.add!(primitive_dry_coupled.output, Speedy.SoilTemperatureOutput())
# initialize coupled simulation
sim_coupled = Speedy.initialize!(primitive_dry_coupled)
# Speedy.timestep!(sim)
# spin up
Speedy.run!(sim_coupled, period=Day(2), output=true)
# run for a few more days
Speedy.run!(sim_coupled, period=Day(5), output=true)

with_theme(fontsize=18) do
    Tair_fig = heatmap(sim_coupled.diagnostic_variables.grid.temp_grid[:,8] .- 273.15, title="", size=(800,400))
    Tsoil_fig = heatmap(RingGrids.Field(interior(land_state.state.temperature)[:,1,end-4], grid), title="", size=(800,400))
    save("plots/speedy_primitive_dry_coupled_tair_R48.png", Tair_fig, px_per_unit=1)
    save("plots/speedy_primitive_dry_coupled_tsoil_R48.png", Tsoil_fig, px_per_unit=1)
    Tair_fig
end

# animate surface layer air temperature
Speedy.animate(sim, variable="temp", coastlines=false, level=spectral_grid.nlayers, output_file="plots/speedy_terrarium_dry_air_temperature.mp4")
Speedy.animate(sim, variable="st", coastlines=false, level=1, output_file="plots/speedy_terrarium_dry_soil_temperature.mp4")

@show land_state.state.air_temperature
@show land_state.state.skin_temperature
@show sim.prognostic_variables.land.soil_temperature

heatmap(sim.diagnostic_variables.grid.temp_grid[:,end])
heatmap(sim.prognostic_variables.land.soil_temperature[:,end])

# Speedy standalone

spectral_grid2 = Speedy.SpectralGrid(trunc=31)
primitive_dry = Speedy.PrimitiveDryModel(spectral_grid2)
sim2 = Speedy.initialize!(primitive_dry)
Speedy.run!(sim2)

using BenchmarkTools

@profview @btime Speedy.timestep!($sim)

default_spectral_grid = Speedy.SpectralGrid(grid.rings)
default_primitive_dry = Speedy.PrimitiveDryModel(default_spectral_grid)
default_sim = Speedy.initialize!(default_primitive_dry)
Speedy.run!(default_sim, period=Day(1))
