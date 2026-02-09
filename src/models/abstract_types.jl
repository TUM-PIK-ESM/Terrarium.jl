# AbstractModel subtypes

# TODO: define general method interfaces (as needed) for all model types

"""
    $TYPEDEF
    
Base type for ground (e.g. soil and rock) models.
"""
abstract type AbstractGroundModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for soil ground models.
"""
abstract type AbstractSoilModel{NF, GR} <: AbstractGroundModel{NF, GR} end

"""
    $TYPEDEF

Base type for land-atmosphere energy exchange models.
"""
abstract type AbstractSurfaceEnergyModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for snow models.
"""
abstract type AbstractSnowModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for vegetation models.
"""
abstract type AbstractVegetationModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    $TYPEDEF

Base type for surface hydrology models.
"""
abstract type AbstractHydrologyModel{NF, GR} <: AbstractModel{NF, GR} end

"""
    AbstractLandModel <: AbstractModel

Base type for full land models which couple together multiple component models.
"""
abstract type AbstractLandModel{NF, GR} <: AbstractModel{NF, GR} end
