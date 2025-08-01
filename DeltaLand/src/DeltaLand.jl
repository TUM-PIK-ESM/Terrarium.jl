module DeltaLand

using DocStringExtensions

import ConstructionBase: getproperties

import DataStructures: OrderedDict

# Oceananigans numerics
# TODO: Raise an issue on Oceananigans.jl about refactoring numerics
# into a separate package.
import Oceananigans
import Oceananigans.Architectures: AbstractArchitecture, CPU
import Oceananigans.BoundaryConditions: fill_halo_regions!
import Oceananigans.Grids: Grids, Periodic, Flat, Bounded
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
export UniformSpacing, ExponentialSpacing, ManualSpacing
include("grids/vertical_discretization.jl")

export ColumnGrid, GlobalRingGrid
include("grids/grids.jl")

# timestepping
include("timesteppers/abstract_timestepper.jl")

export ForwardEuler
include("timesteppers/forward_euler.jl")

export get_grid, get_time_stepping, initialize, update_state!, compute_tendencies!, timestep!
include("models/abstract_model.jl")

# physical processes
include("processes/processes.jl")

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
