module Terrarium

using DocStringExtensions

import Base: @propagate_inbounds

import ConstructionBase: ConstructionBase, getproperties, setproperties

import DataStructures: OrderedDict

import Dates: Dates, TimeType, Period, Year, Month, Day, Hour, Minute, Second

import DomainSets: RealLine, HalfLine, PositiveRealLine, UnitInterval, AbstractInterval

import Flatten: flatten, flattenable, reconstruct

import Interpolations

# Oceananigans numerics
import Oceananigans
import Oceananigans.AbstractOperations: Average, Integral
import Oceananigans.Advection: AbstractAdvectionScheme, UpwindBiased
import Oceananigans.Architectures: AbstractArchitecture, CPU, GPU, architecture, on_architecture, array_type
import Oceananigans.Fields: Field, FunctionField, AbstractField, Center, Face, set!, interior, xnodes, ynodes, znodes, location
import Oceananigans.Grids as OceananigansGrids
import Oceananigans.Grids: Periodic, Flat, Bounded
import Oceananigans.Operators: ∂zᵃᵃᶜ, ∂zᵃᵃᶠ, ℑzᵃᵃᶠ, Δzᵃᵃᶜ
import Oceananigans.OutputReaders: FieldTimeSeries
import Oceananigans.Simulations: Simulation, run!, timestepper
import Oceananigans.TimeSteppers: Clock, update_state!, time_step!, tick!, reset!
import Oceananigans.Units: Time
import Oceananigans.Utils: launch!
# Boundary conditions
import Oceananigans.BoundaryConditions: BoundaryCondition, DefaultBoundaryCondition, FieldBoundaryConditions,
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

"""
Alias for numeric `Quantity` with type `NF` and units `U`.
"""
const LengthQuantity{NF, U} = Quantity{NF, 𝐋, U} where {NF, U<:Units}

"""
Alias for Oceananigans `AbstractBoundaryConditionClassification`.
"""
const BCType = AbstractBoundaryConditionClassification

# Re-export selected types and methods from Oceananigans
export Simulation, Field, FieldTimeSeries, CPU, GPU, Clock, Center, Face
export Value, Flux, Gradient, ValueBoundaryCondition, GradientBoundaryCondition, FluxBoundaryCondition, NoFluxBoundaryCondition
export run!, time_step!, set!, interior, architecture, on_architecture, xnodes, ynodes, znodes, location

# Re-export common Dates types
export Year, Month, Day, Hour, Minute, Second

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

export InputSource, InputSources, FieldInputSource, FieldTimeSeriesInputSource
export update_inputs!
include("input_output/input_sources.jl")

# timestepping
export timestep!, default_dt, is_adaptive
include("timesteppers/abstract_timestepper.jl")

# process/model interface
export get_grid, get_initializer, variables, compute_auxiliary!, compute_tendencies!
include("abstract_model.jl")

# state variables
export StateVariables, get_fields
include("state_variables.jl")

# default initializers
export FieldInitializers, DefaultInitializer
include("initializers.jl")

export FieldBC, FieldBCs, boundary_conditions
include("boundary_conditions.jl")

# physical processes
include("processes/processes.jl")

# concrete model implementations
include("models/models.jl")

# timestepper implementations
export ForwardEuler
include("timesteppers/forward_euler.jl")
export Heun
include("timesteppers/heun.jl")

# model integrator/simulation types and methods
export ModelIntegrator, initialize, current_time, iteration
include("timesteppers/model_integrator.jl")

end # module Terrarium
