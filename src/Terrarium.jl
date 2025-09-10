module Terrarium

using DocStringExtensions

import ConstructionBase: getproperties, setproperties

import DataStructures: OrderedDict
import Dates: Period, Second

# Oceananigans numerics
# TODO: Raise an issue on Oceananigans.jl about refactoring numerics
# into a separate package.
import Oceananigans: Field, Center, Face, set!
import Oceananigans.Advection: AbstractAdvectionScheme, UpwindBiased
import Oceananigans.Architectures: AbstractArchitecture, CPU, GPU, architecture, on_architecture
import Oceananigans.BoundaryConditions: BoundaryConditions, fill_halo_regions!
import Oceananigans.Grids: Grids, Periodic, Flat, Bounded
import Oceananigans.Operators: ∂zᵃᵃᶜ, ∂zᵃᵃᶠ, ℑzᵃᵃᶠ, Δzᵃᵃᶜ
import Oceananigans.TimeSteppers: Clock, tick_time!, reset!
import Oceananigans.Utils: launch!

# Adapt and KernelAbstractions for GPU parallelization
import Adapt: Adapt, adapt
import KernelAbstractions: @kernel, @index

# Freeze curves for soil energy balance
import FreezeCurves

import RingGrids

# Re-export important types and methods
export CPU, GPU, Clock
export adapt, set!

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

# default initializers
export FieldInitializers
include("models/initializers.jl")

export FieldBoundaryConditions
include("models/boundary_conditions.jl")

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
export Simulation, initialize, run!, current_time
include("simulation.jl")

end # module Terrarium
