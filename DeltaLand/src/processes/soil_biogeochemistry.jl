# Temporary stub implementation of BGC
Base.@kwdef struct ConstantSoilCarbonDenisty{NF} <: AbstractSoilBiogeochemistry
    "Soil organic carbon density [kg/m^3]"
    ρ_soc::NF = 65.0

    "Pure organic matter density [kg/m^3]"
    ρ_org::NF = 1300.0
end

@inline update_state(i, j, k, state, model, strat::ConstantSoilCarbonDenisty) = nothing

@inline compute_tendencies(i, j, k, state, model, strat::ConstantSoilCarbonDenisty) = nothing
