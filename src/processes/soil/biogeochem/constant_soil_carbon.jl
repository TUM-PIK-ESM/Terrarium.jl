"""
    $TYPEDEF

Naive implementation of soil biogeochemistry that just assumes there to be a constant
organic content in all soil layers.

Properties:
$TYPEDFIELDS
"""
@parameterized @kwdef struct ConstantSoilCarbonDensity{NF} <: AbstractSoilBiogeochemistry{NF}
    "Soil organic carbon density"
    @param ρ_soc::NF = 0.0 (units = u"kg/m^3", bounds = Positive)

    "Pure organic matter density"
    @param ρ_org::NF = 1300.0 (units = u"kg/m^3", bounds = Positive)
end

ConstantSoilCarbonDensity(::Type{NF}; kwargs...) where {NF} = ConstantSoilCarbonDensity{NF}(; kwargs...)

density_pure_soc(bgc::ConstantSoilCarbonDensity) = bgc.ρ_org

variables(::ConstantSoilCarbonDensity) = ()

"""
    $SIGNATURES

Calculate the organic solid fraction based on the prescribed SOC and natural porosity/density of
the organic material.
"""
@propagate_inbounds density_soc(i, j, k, grid, fields, bgc::ConstantSoilCarbonDensity) = bgc.ρ_soc

@inline compute_auxiliary!(state, grid, bgc::ConstantSoilCarbonDensity, args...) = nothing

@inline compute_tendencies!(state, grid, bgc::ConstantSoilCarbonDensity, args...) = nothing
