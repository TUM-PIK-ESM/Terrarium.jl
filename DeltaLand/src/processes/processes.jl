include("abstract_types.jl")

export PhysicalConstants
include("physical_constants.jl")

# Soil

export SoilThermalConductivities, SoilHeatCapacities, SoilThermalProperties, InverseQuadratic
include("soil/soil_thermal_properties.jl")

export SURFEXHydraulics
include("soil/soil_hydraulic_properties.jl")

export HomogeneousSoil
include("soil/soil_stratigraphy.jl")

export ImmobileSoilWater, SoilHydraulicProperties
include("soil/soil_hydrology.jl")

export SoilEnergyBalance
include("soil/soil_energy.jl")

export ConstantSoilCarbonDenisty
include("soil/soil_biogeochemistry.jl")

# Vegetation

include("vegetation/photosynthesis.jl")

include("vegetation/carbon_dynamics.jl")
