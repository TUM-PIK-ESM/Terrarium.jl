module DeltaLand

using DocStringExtensions

import ConstructionBase: getproperties

import DataStructures: OrderedDict

# Oceananigans numerics
using Oceananigans
import Oceananigans.Architectures: AbstractArchitecture
import Oceananigans.BoundaryConditions: fill_halo_regions!
import Oceananigans.Operators: ∂zᵃᵃᶜ, ∂zᵃᵃᶠ, Δzᵃᵃᶜ
import Oceananigans.TimeSteppers: Clock, tick_time!, reset!
import Oceananigans.Utils: launch!

# KernelAbstractions for GPU parallelization
import KernelAbstractions: @kernel, @index

# Freeze curves for soil energy balance
import FreezeCurves

# temporary dependency on SpeedyWeather until RingGrids is registered
import SpeedyWeather: RingGrids

# temporary dependency on CryoGrid for soil types and SEB
import CryoGrid: SoilTexture, SurfaceEnergyBalance

# internal utilities
include("utils.jl")

# grids
export GlobalRingGrid, UniformSpacing, ExponentialSpacing, ManualSpacing, ModelInitializer
include("grids/vertical_discretization.jl")
include("grids/grids.jl")
include("grids/initializers.jl")

# timestepping
include("timesteppers/abstract_timestepper.jl")

export ForwardEuler
include("timesteppers/forward_euler.jl")

# abstract types
include("processes/abstract_types.jl")

export get_grid, get_time_stepping, initialize, update_state!, compute_tendencies!, timestep!
include("models/abstract_model.jl")

# physical processes

export PhysicalConstants
include("processes/physical_constants.jl")

export ImmobileSoilWater, SoilHydraulicProperties
include("processes/soil_hydrology.jl")

export SoilThermalConductivities, SoilHeatCapacities
include("processes/soil_thermal_properties.jl")

export SoilEnergyBalance, SoilThermalProperties
include("processes/soil_energy.jl")

export HomogeneousSoil
include("processes/soil_stratigraphy.jl")

export ConstantSoilCarbonDenisty
include("processes/soil_biogeochemistry.jl")

# state variables
include("state_variables.jl")

# simulation types
export Simulation
include("simulation.jl")

# concrete model implementations
export SoilModel
include("models/soil_model.jl")

export LandModel
include("models/land_model.jl")

end # module DeltaLand
