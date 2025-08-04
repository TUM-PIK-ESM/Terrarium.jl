module DeltaLand

using DocStringExtensions

import ConstructionBase: getproperties

import DataStructures: OrderedDict

# Oceananigans numerics
# TODO: Raise an issue on Oceananigans.jl about refactoring numerics
# into a separate package.
import Oceananigans: Field, Center, Face
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

export PrognosticVariable, AuxiliaryVariable
include("abstract_variables.jl")

# grids
export UniformSpacing, ExponentialSpacing, PrescribedSpacing
include("grids/vertical_discretization.jl")

export ColumnGrid, GlobalRingGrid, get_field_grid
include("grids/grids.jl")

# timestepping
export timestep!, get_dt, is_adaptive
include("timesteppers/abstract_timestepper.jl")

# model interface
export get_grid, get_time_stepping, get_boundary_conditions, variables, compute_auxiliary!, compute_tendencies!
include("models/abstract_model.jl")

# timestepper implementations
export ForwardEuler
include("timesteppers/forward_euler.jl")

# default boundary conditions
export FieldBoundaryConditions, PrescribedFluxes
include("models/boundary_conditions.jl")

# default initializers
export FieldInitializers
include("models/initializers.jl")

# state variables
export StateVariables
include("state_variables.jl")

# physical processes
include("processes/processes.jl")

# concrete model implementations
export SoilModel
include("models/soil_model.jl")

export VegetationModel
include("models/vegetation_model.jl")

export LandModel
include("models/land_model.jl")

# simulation types
export Simulation, initialize
include("simulation.jl")

end # module DeltaLand
