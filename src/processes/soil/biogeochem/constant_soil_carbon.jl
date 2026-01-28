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
end

ConstantSoilCarbonDensity(::Type{NF}; kwargs...) where {NF} = ConstantSoilCarbonDensity{NF}(; kwargs...)

variables(::ConstantSoilCarbonDensity) = ()

"""
    $SIGNATURES

Calculate the organic solid fraction based on the prescribed SOC and natural porosity/density of
the organic material.
"""
@propagate_inbounds organic_fraction(bgc::ConstantSoilCarbonDensity) = bgc.ρ_soc / ((1 - bgc.por_org)*bgc.ρ_org)
@propagate_inbounds organic_fraction(i, j, k, grid, state, bgc::ConstantSoilCarbonDensity) = organic_fraction(bgc)

@inline initialize!(state, model, bgc::ConstantSoilCarbonDensity) = nothing

@inline compute_auxiliary!(state, grid, bgc::ConstantSoilCarbonDensity) = nothing

@inline compute_tendencies!(state, grid, bgc::ConstantSoilCarbonDensity) = nothing
