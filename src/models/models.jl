# Abstract types

# TODO: define general method interfaces (as needed) for all model types

"""
    $TYPEDEF
    
Base type for ground (e.g. soil and rock) models.
"""
abstract type AbstractGroundModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for soil ground models.
"""
abstract type AbstractSoilModel{NF, GR, TS} <: AbstractGroundModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for snow models.
"""
abstract type AbstractSnowModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for vegetation/carbon models.
"""
abstract type AbstractVegetationModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for surface energy balance models.
"""
abstract type AbstractSurfaceEnergyBalanceModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    $TYPEDEF

Base type for surface hydrology models.
"""
abstract type AbstractHydrologyModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

"""
    AbstractLandModel <: AbstractModel

Base type for full land models which couple together multiple component models.
"""
abstract type AbstractLandModel{NF, GR, TS} <: AbstractModel{NF, GR, TS} end

# Soil

export SoilModel
include("soil/soil_model.jl")

export SoilBoundaryConditions, SoilBC, GroundHeatFlux, GeothermalHeatFlux, FreeDrainage, ImpermeableBoundary
include("soil/soil_model_bcs.jl")

include("soil/soil_model_init.jl")

# Vegetation

export VegetationModel
include("vegetation/vegetation_model.jl")

# Coupled Land

export LandModel
include("land/land_model.jl")
