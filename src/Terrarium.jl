module Terrarium

using DocStringExtensions

using Adapt: Adapt, adapt, @adapt_structure

using Base: @propagate_inbounds

using ConstructionBase: ConstructionBase, getproperties, setproperties

using DataStructures: OrderedDict

using Dates: Dates, TimeType, Period, Year, Month, Day, Hour, Minute, Second

using DomainSets: RealLine, HalfLine, PositiveRealLine, UnitInterval, AbstractInterval

using Flatten: flatten, flattenable, reconstruct

using Interpolations

using KernelAbstractions: @kernel, @index

# Oceananigans numerics
using Oceananigans.AbstractOperations: Average, Integral
using Oceananigans.Architectures: Architectures, AbstractArchitecture, CPU, GPU, architecture, on_architecture, array_type
using Oceananigans.Fields: Field, FunctionField, AbstractField, Center, Face, set!, compute!, interior, xnodes, ynodes, znodes, zspacings, location
using Oceananigans.Forcings: Forcings, Forcing, ContinuousForcing, DiscreteForcing
using Oceananigans.Grids: Periodic, Flat, Bounded
using Oceananigans.Operators: ∂zᵃᵃᶜ, ∂zᵃᵃᶠ, ℑzᵃᵃᶠ, Δzᵃᵃᶜ
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.Simulations: Simulation, run!, timestepper
using Oceananigans.TimeSteppers: Clock, update_state!, time_step!, tick!, reset!
using Oceananigans.Units: Time
using Oceananigans.Utils: launch!

# Boundary conditions
using Oceananigans.BoundaryConditions: BoundaryConditions, BoundaryCondition, DefaultBoundaryCondition, FieldBoundaryConditions,
                                       ValueBoundaryCondition, FluxBoundaryCondition, GradientBoundaryCondition, NoFluxBoundaryCondition,
                                       ContinuousBoundaryFunction, DiscreteBoundaryFunction,
                                       AbstractBoundaryConditionClassification, Value, Flux, Gradient, # BC type classifications
                                       fill_halo_regions!, regularize_field_boundary_conditions, getbc, compute_z_bcs!

# Freeze curves for soil energy balance
using FreezeCurves: FreezeCurves, FreezeCurve, SFCC, SWRC, FreeWater

# Units (for testing and UI)
# Unit dimensions for length (𝐋), mass (𝐌), and time (𝐓)
using Unitful: 𝐋, 𝐌, 𝐓
using Unitful: Units, Quantity, AbstractQuantity, NoUnits
using Unitful: @u_str, uconvert, ustrip, upreferred

# Explicit imports
import Oceananigans
import RingGrids

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
export run!, time_step!, set!, compute!, interior, architecture, on_architecture, xnodes, ynodes, znodes, zspacings, location

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
