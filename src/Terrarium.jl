module Terrarium

using DocStringExtensions

import ConstructionBase: getproperties, setproperties

import DataStructures: OrderedDict

import Dates: Period, Year, Month, Day, Second

import Flatten

import Interpolations

import NamedTupleTools: merge_recursive

# Oceananigans numerics
# TODO: Raise an issue on Oceananigans.jl about refactoring numerics
# into a separate package.
import Oceananigans: Field, AbstractField, Center, Face, set!, interior, xnodes, ynodes, znodes
import Oceananigans.Advection: AbstractAdvectionScheme, UpwindBiased
import Oceananigans.Architectures: AbstractArchitecture, CPU, GPU, architecture, on_architecture
import Oceananigans.Grids: Grids, Periodic, Flat, Bounded
import Oceananigans.Operators: ∂zᵃᵃᶜ, ∂zᵃᵃᶠ, ℑzᵃᵃᶠ, Δzᵃᵃᶜ
import Oceananigans.TimeSteppers: Clock, tick_time!, reset!
import Oceananigans.Units: Time
import Oceananigans.Utils: launch!
# Boundary conditions
import Oceananigans.BoundaryConditions: FieldBoundaryConditions, BoundaryCondition, DefaultBoundaryCondition,
                                        ValueBoundaryCondition, FluxBoundaryCondition, NoFluxBoundaryCondition,
                                        fill_halo_regions!, regularize_field_boundary_conditions

# Adapt and KernelAbstractions for GPU parallelization
import Adapt: Adapt, adapt
import KernelAbstractions: @kernel, @index

# Freeze curves for soil energy balance
import FreezeCurves

# temporary dependency on SpeedyWeather until RingGrids is registered
import SpeedyWeather: RingGrids

# Units (for testing and UI)
# Unit dimensions for length (𝐋), mass (𝐌), and time (𝐓)
import Unitful: 𝐋, 𝐌, 𝐓
import Unitful: Units, Quantity, AbstractQuantity, NoUnits
import Unitful: @u_str, uconvert, ustrip, upreferred

const LengthQuantity{NF, U} = Quantity{NF, 𝐋, U} where {NF, U<:Units}

# Re-export selected types and methods from Oceananigans
export CPU, GPU, Clock, Center, Face, ValueBoundaryCondition, FluxBoundaryCondition, NoFluxBoundaryCondition
export set!, interior, architecture, on_architecture, xnodes, ynodes, znodes

# Re-export Dates types
export Year, Month, Day, Second

# Re-export adapt
export adapt

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
include("abstract_model.jl")

# default initializers
export Initializers, DefaultInitializer
include("initializers.jl")

export VerticalBoundaryConditions, DefaultBoundaryConditions, PrescribedFlux
include("boundary_conditions.jl")

# timestepper implementations
export ForwardEuler
include("timesteppers/forward_euler.jl")

# state variables
export StateVariables
include("state_variables.jl")

# physical processes
include("processes/processes.jl")

# concrete model implementations
include("models/models.jl")

# simulation types
export Simulation, initialize, run!, current_time
include("simulation.jl")

end # module Terrarium
