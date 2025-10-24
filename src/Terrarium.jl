module Terrarium

using DocStringExtensions

import ConstructionBase: getproperties, setproperties

import DataStructures: OrderedDict

import Dates: Dates, TimeType, Period, Year, Month, Day, Second

import DomainSets: RealLine, HalfLine, PositiveRealLine, UnitInterval, AbstractInterval

import Flatten: flatten, flattenable, reconstruct

import Interpolations

# Oceananigans numerics
# TODO: Raise an issue on Oceananigans.jl about refactoring numerics into a separate package.
import Oceananigans: Oceananigans, Field, AbstractField, Center, Face, set!, interior, xnodes, ynodes, znodes, location
import Oceananigans.Advection: AbstractAdvectionScheme, UpwindBiased
import Oceananigans.Architectures: AbstractArchitecture, CPU, GPU, architecture, on_architecture, array_type
import Oceananigans.Grids as OceananigansGrids
import Oceananigans.Grids: Periodic, Flat, Bounded
import Oceananigans.Operators: ∂zᵃᵃᶜ, ∂zᵃᵃᶠ, ℑzᵃᵃᶠ, Δzᵃᵃᶜ
import Oceananigans.OutputReaders: FieldTimeSeries
import Oceananigans.Simulations: Simulation, run!, timestepper
import Oceananigans.TimeSteppers: Clock, update_state!, time_step!, tick!, reset!
import Oceananigans.Units: Time
import Oceananigans.Utils: launch!
# Boundary conditions
import Oceananigans.BoundaryConditions: FieldBoundaryConditions, BoundaryCondition, DefaultBoundaryCondition,
                                        ValueBoundaryCondition, FluxBoundaryCondition, GradientBoundaryCondition, NoFluxBoundaryCondition,
                                        ContinuousBoundaryFunction, DiscreteBoundaryFunction,
                                        AbstractBoundaryConditionClassification, Value, Flux, Gradient, # BC type classifications
                                        fill_halo_regions!, regularize_field_boundary_conditions, getbc, compute_z_bcs!

# Adapt and KernelAbstractions for GPU parallelization
import Adapt: Adapt, adapt, @adapt_structure
import KernelAbstractions: @kernel, @index

# Freeze curves for soil energy balance
import FreezeCurves: FreezeCurves, FreezeCurve, SFCC, SWRC, FreeWater

import RingGrids

# Units (for testing and UI)
# Unit dimensions for length (𝐋), mass (𝐌), and time (𝐓)
import Unitful: 𝐋, 𝐌, 𝐓
import Unitful: Units, Quantity, AbstractQuantity, NoUnits
import Unitful: @u_str, uconvert, ustrip, upreferred

const LengthQuantity{NF, U} = Quantity{NF, 𝐋, U} where {NF, U<:Units}

const BCType = AbstractBoundaryConditionClassification

# Re-export selected types and methods from Oceananigans
export Simulation, Field, FieldTimeSeries, CPU, GPU, Clock, Center, Face
export Value, Flux, Gradient, ValueBoundaryCondition, GradientBoundaryCondition, FluxBoundaryCondition, NoFluxBoundaryCondition
export run!, time_step!, set!, interior, architecture, on_architecture, xnodes, ynodes, znodes, location

# Re-export common Dates types
export Year, Month, Day, Second

# Re-export unit types
export @u_str, uconvert, ustrip

# Re-export adapt
export adapt

# internal utilities
include("utils.jl")

export PrognosticVariable, AuxiliaryVariable, InputVariable, Input, XY, XYZ
include("abstract_variables.jl")

# grids
export UniformSpacing, ExponentialSpacing, PrescribedSpacing
include("grids/vertical_discretization.jl")

export ColumnGrid, ColumnRingGrid, get_field_grid
include("grids/grids.jl")

export InputFields, InputProvider, InputSource
export update_inputs!, get_input_fields, get_input_field
include("inputs/inputs.jl")

# timestepping
export timestep!, default_dt, is_adaptive
include("timesteppers/abstract_timestepper.jl")

# model interface
export get_grid, get_time_stepping, get_boundary_conditions, variables, compute_auxiliary!, compute_tendencies!
include("abstract_model.jl")

# state variables
export StateVariables, get_fields
include("state_variables.jl")

# default initializers
export FieldInitializers, DefaultInitializer
include("initializers.jl")

export ColumnBoundaryConditions, ColumnBCs, DefaultBoundaryConditions, PrescribedFlux, PrescribedValue, PrescribedGradient
include("boundary_conditions.jl")

# timestepper implementations
export ForwardEuler
include("timesteppers/forward_euler.jl")

# physical processes
include("processes/processes.jl")

# concrete model implementations
include("models/models.jl")

# model simulation types and methods
export ModelState, initialize, current_time, iteration
include("model_state.jl")

end # module Terrarium
