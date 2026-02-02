# Vegetation process types

abstract type AbstractPhotosynthesis{NF} <: AbstractProcess{NF} end

abstract type AbstractStomatalConductance{NF} <: AbstractProcess{NF} end

abstract type AbstractAutotrophicRespiration{NF} <: AbstractProcess{NF} end

abstract type AbstractPlantAvailableWater{NF} <: AbstractProcess{NF} end

abstract type AbstractVegetationDynamics{NF} <: AbstractProcess{NF} end

abstract type AbstractPhenology{NF} <: AbstractProcess{NF} end

abstract type AbstractVegetationCarbonDynamics{NF} <: AbstractProcess{NF} end

# Parameterizations

abstract type AbstractRootDistribution{NF} end
