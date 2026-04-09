# Surface hydrology process types

"""
Base type for coupled surface hydrology processes.
"""
abstract type AbstractSurfaceHydrology{NF} <: AbstractCoupledProcesses{NF} end

get_evapotranspiration(hydrology::AbstractSurfaceHydrology) = hydrology.evapotranspiration

include("canopy_interception/abstract_types.jl")
include("evapotranspiration/abstract_types.jl")
include("runoff/abstract_types.jl")
