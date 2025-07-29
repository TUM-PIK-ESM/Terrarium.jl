module DeltaLand

using DocStringExtensions

import DataStructures: OrderedDict

# Oceananigans numerics
import Oceananigans
import Oceananigans.BoundaryConditions: fill_halo_regions!
import Oceananigans.Operators
import Oceananigans.TimeSteppers: Clock
import Oceananigans.Utils: launch!

# KernelAbstractions for GPU parallelization
import KernelAbstractions: @kernel, @index

# temporary dependency on SpeedyWeather until RingGrids is registered
import SpeedyWeather: RingGrids

# temporary dependency on CryoGrid for soil types and SEB
import CryoGrid: SoilTexture, SurfaceEnergyBalance

# grids
export GlobalRingGrid, UniformSpacing, ExponentialSpacing, ManualSpacing
include("grids/grids.jl")

# timestepping
export Euler
include("timesteppers/abstract_timestepper.jl")
include("timesteppers/forward_euler.jl")

# abstract types
include("processes/abstract_types.jl")

export get_grid, get_time_stepping, initialize, update_state!, compute_tendencies!, timestep!
include("models/abstract_model.jl")

# physical processes

export PhyscialConstants
include("processes/physical_constants.jl")

export ImmobileSoilWater, SoilHydraulicProperties
include("processes/soil_hydrology.jl")

export SoilEnergyBalance, SoilThermalProperties
include("processes/soil_energy.jl")

export HomogeneousSoil
include("processes/soil_stratigraphy.jl")

export ConstantSoilCarbonDenisty
include("processes/soil_biogeochemistry.jl")

# simulation types
export Simulation
include("simulation.jl")

# concrete model implementations
export SoilModel
include("models/soil_model.jl")

export LandSurfaceModel
include("models/land_surface_model.jl")

export LandModel
include("models/land_model.jl")

end # module DeltaLand
