"""
    $TYPEDEF

Naive implementation of soil biogeochemistry that just assumes there to be a constant
organic content in all soil layers.

Properties:
$TYPEDFIELDS
"""
Base.@kwdef struct ConstantSoilCarbonDensity{NF} <: AbstractSoilBiogeochemistry{NF}
    "Soil organic carbon density [kg/m^3]"
    ρ_soc::NF = 0.0

    "Pure organic matter density [kg/m^3]"
    ρ_org::NF = 1300.0

    "Natural porosity of organic material"
    por_org::NF = 0.9
end

ConstantSoilCarbonDensity(::Type{NF}; kwargs...) where {NF} = ConstantSoilCarbonDensity{NF}(; kwargs...)

variables(::ConstantSoilCarbonDensity) = ()

"""
    organic_porosity(i, j, k, state, bgc::ConstantSoilCarbonDensity)

Get the prescribed natural porosity of organic soil.
"""
@inline organic_porosity(i, j, k, state, bgc::ConstantSoilCarbonDensity) = bgc.por_org

"""
    organic_fraction(i, j, k, state, bgc::ConstantSoilCarbonDensity)

Calculate the organic solid fraction at the given indices based on the prescribed SOC and natural porosity/density of
the organic material.
"""
@inline organic_fraction(i, j, k, state, bgc::ConstantSoilCarbonDensity) = bgc.ρ_soc / ((1 - bgc.por_org) * bgc.ρ_org)

@inline initialize!(state, model, bgc::ConstantSoilCarbonDensity) = nothing

@inline compute_auxiliary!(state, model, bgc::ConstantSoilCarbonDensity) = nothing

@inline compute_tendencies!(state, model, bgc::ConstantSoilCarbonDensity) = nothing
