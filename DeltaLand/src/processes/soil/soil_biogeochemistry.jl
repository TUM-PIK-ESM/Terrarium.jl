# Temporary stub implementation of BGC
Base.@kwdef struct ConstantSoilCarbonDenisty{NF} <: AbstractSoilBiogeochemistry
    "Soil organic carbon density [kg/m^3]"
    ρ_soc::NF = 65.0

    "Pure organic matter density [kg/m^3]"
    ρ_org::NF = 1300.0

    "Natural porosity of organic material"
    por_org::NF = 0.90
end

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

@inline organic_fraction(idx, state, bgc::ConstantSoilCarbonDenisty) = organic_fraction(bgc)

@inline compute_auxiliary!(state, model, bgc::ConstantSoilCarbonDenisty) = nothing

@inline compute_tendencies!(state, model, bgc::ConstantSoilCarbonDenisty) = nothing
