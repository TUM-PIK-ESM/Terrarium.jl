# Note: Should sand, silt, and clay have separate thermal properties?
@kwdef struct SoilThermalConductivities{NF}
    water::NF = 0.57 # thermal conductivity of water [W/m/K] (Hillel 1982)
    ice::NF = 2.2 # thermal conductivity of ice [W/m/K] Hillel (1982)
    air::NF = 0.025 # thermal conductivity of air [W/m/K] Hillel (1982)
    mineral::NF = 3.8 # thermal conductivity of mineral soil constituents [W/m/K] Hillel (1982)
    organic::NF = 0.25 # thermal conductivity of organic soil constituents [W/m/K] Hillel (1982)
end

"""
Base type for bulk weighting/mixing schemes that calculate weighted mixture of material properties
such as conductivities or densities.
"""
abstract type AbstractBulkWeightingScheme end

@kwdef struct SoilHeatCapacities{NF}
    water::NF = 4.2e6 # volumetric heat capacity of water [J/m^3]
    ice::NF = 1.9e6 # volumetric heat capacity of ice [J/m^3]
    air::NF = 0.00125e6 # volumetric heat capacity of air [J/m^3]
    mineral::NF = 2.0e6 # volumetric heat capacity of mineral soil [J/m^3]
    organic::NF = 2.5e6 # volumetric heat capacity of organic soil [J/m^3]
end

# TODO: In principle, these types could change for different soil parameterizations.
# This is something we should ideally allow for.
@kwdef struct SoilThermalProperties{NF,CondWeighting}
    "Thermal conductivities for all constituents"
    cond::SoilThermalConductivities{NF} = SoilThermalConductivities()

    "Thermal conductivity mixing scheme"
    cond_bulk::CondWeighting = InverseQuadratic()

    "Thermal conductivities for all constituents"
    heatcap::SoilHeatCapacities{NF} = SoilHeatCapacities()
end

"""
The inverse quadratic (or "quadratic parallel") bulk thermal conductivity formula (Cosenza et al. 2003):

```math
k = [\\sum_{i=1}^N θᵢ\\sqrt{kᵢ}]^2
```

Cosenza, P., Guérin, R., and Tabbagh, A.: Relationship between thermal
conductivity and water content of soils using numerical modelling,
European Journal of Soil Science, 54, 581–588,
https://doi.org/10.1046/j.1365-2389.2003.00539.x, 2003.
"""
struct InverseQuadratic <: AbstractBulkWeightingScheme end

(f::InverseQuadratic)(x::Real, weight::Real) = sqrt(x)*weight
# we use fastmap here so that the ordering of named tuples can be arbitrary
(f::InverseQuadratic)(xs, weights) = sum(fastmap(f, xs, weights))^2
