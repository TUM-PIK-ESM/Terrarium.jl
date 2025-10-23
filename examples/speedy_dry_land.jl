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
        return new{typeof(model)}(spectral_grid, model)
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
    Tair = diagn.grid.temp_grid[:,end] .- NF(273.15)
    Tair_land = get_input_field(land.model.inputs.fields, :air_temperature, XY())
    set!(Tair_land, Tair)
    # do timestep
    Terrarium.timestep!(land.model, Dates.second(progn.clock.Δt))
    # extract surface temperature and copy to speedy prognostic variable
    Tsoil = land.model.state.temperature
    Tsurf = Tsoil[1:Tsoil.grid.Nx, 1, Tsoil.grid.Nz:Tsoil.grid.Nz]
    progn.land.soil_temperature .= Tsurf .+ NF(273.15)
end

ring_grid = RingGrids.FullGaussianGrid(24)
grid = ColumnRingGrid(CPU(), Float32, ExponentialSpacing(N=30), ring_grid)
model = SoilModel(grid)
# Initial conditions
initializer = FieldInitializers(
    # steady-ish state initial condition for soil temperature
    temperature = (x,z) -> 10 + 0.02*z,
    # fully saturated soil
    saturation_water_ice = 1.0,
)
# Periodic surface temperature with annual cycle
T_ub = PrescribedValue(:temperature, Input(:air_temperature, units=u"°C"))
boundary_conditions = SoilBoundaryConditions(eltype(grid), top=T_ub)
model = SoilModel(grid; initializer, boundary_conditions)
state = initialize(model)

# Set up Speedy land and atmosphere model
land = TerrariumDryLand(state)
spectral_grid = land.spectral_grid #Speedy.SpectralGrid(grid.rings)
land_sea_mask = Speedy.RockyPlanetMask(spectral_grid)
output = Speedy.NetCDFOutput(spectral_grid, Speedy.PrimitiveDryModel, path="outputs/")
primitive_dry = Speedy.PrimitiveDryModel(spectral_grid; land, land_sea_mask, output)
Speedy.add!(primitive_dry.output, Speedy.SoilTemperatureOutput())
sim = Speedy.initialize!(primitive_dry)
Speedy.run!(sim, period=Day(30), output=true)
# animate surface layer air temperature
Speedy.animate(sim, variable="temp", coastlines=false, level=spectral_grid.nlayers, output_file="plots/speedy_terrarium_dry_air_temperature.mp4")
Speedy.animate(sim, variable="st", coastlines=false, level=1, output_file="plots/speedy_terrarium_dry_soil_temperature.mp4")

@show state.state.air_temperature
@show state.state.temperature
@show sim.prognostic_variables.land.soil_temperature

heatmap(sim.diagnostic_variables.grid.temp_grid[:,end])
heatmap(sim.prognostic_variables.land.soil_temperature[:,end])

@profview @btime Speedy.timestep!(sim.prognostic_variables, sim.diagnostic_variables, land, sim.model)