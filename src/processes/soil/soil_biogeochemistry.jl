"""
    $TYPEDEF

Naive implementation of soil biogeochemistry that just assumes there to be a constant
organic content in all soil layers.

Properties:
$TYPEDFIELDS
"""
Base.@kwdef struct ConstantSoilCarbonDenisty{NF} <: AbstractSoilBiogeochemistry{NF}
    "Soil organic carbon density [kg/m^3]"
    ρ_soc::NF = 0.0

    "Pure organic matter density [kg/m^3]"
    ρ_org::NF = 1300.0

    "Natural porosity of organic material"
    por_org::NF = 0.90
end

ConstantSoilCarbonDenisty(::Type{NF}; kwargs...) where {NF} = ConstantSoilCarbonDenisty{NF}(; kwargs...)

variables(::ConstantSoilCarbonDenisty) = ()

"""
    organic_porosity(bgc::ConstantSoilCarbonDenisty)

Get the prescribed natural porosity of organic soil.
"""
organic_porosity(bgc::ConstantSoilCarbonDenisty) = bgc.por_org

"""
    organic_fraction(bgc::ConstantSoilCarbonDenisty)

Calculate the organic solid fraction based on the prescribed SOC and natural porosity/density of
the organic material.
"""
organic_fraction(bgc::ConstantSoilCarbonDenisty) = bgc.ρ_soc / ((1 - bgc.por_org)*bgc.ρ_org)

"""
    organic_fraction(idx, state, bgc::ConstantSoilCarbonDenisty)

Calculate the soil organic carbon fraction at the given grid index. For `ConstantSoilCarbonDenisty`,
this simply returns a constant parameter. Implementations of soil carbon dynamics would want to compute
this based on the prognostic state of the carbon content/pool stored in the soil.
"""
@inline organic_fraction(idx, state, bgc::ConstantSoilCarbonDenisty) = organic_fraction(bgc)

@inline compute_auxiliary!(state, model, bgc::ConstantSoilCarbonDenisty) = nothing

@inline compute_tendencies!(state, model, bgc::ConstantSoilCarbonDenisty) = nothing
