@kwdef struct PaladynVegetationDynamics{NF} <: AbstractVegetationDynamics
    "Vegetation seed fraction"
    ν_seed::NF = 0.001

    "Minimum vegetation disturbance rate [1/year]"
    γv_min::NF = 0.002 
end

variables(veg_dynamics::PaladynVegetationDynamics) = (
    prognostic(:ν, XY()), # PFT fractional area coverage
    auxiliary(:C_veg, XY()), # Vegetation carbon pool (kgC/m²)
    auxiliary(:LAI_b, XY()), # Balanced Leaf Area Index (LAI_b)
)

@inline function compute_auxiliary!(idx, state, model, veg_dynamics::PaladynVegetationDynamics)
    # Nothing needed here for now
    return nothing
end

@inline function compute_tendencies!(idx, state, model::AbstractVegetationModel, veg_dynamics::PaladynVegetationDynamics{NF}) where NF
    i, j = idx
    
    # Compute λ_NPP
    λ_NPP = compute_λ_NPP(idx, state, model, model.carbon_dynamics)

    # Compute the disturbance rate
    # TODO add Paladyn implemetation
    # Placeholder for now = min disturbance rate
    γv = veg_dynamics.γv_min

    # Compute the vegetation fraction tendency
    ν_star = max(state.ν[i, j], veg_dynamics.ν_seed) 
    state.ν_tendency[i, j] = (λ_NPP * state.NPP[i, j] / state.C_veg[i, j]) * 
                             ν_star * (NF(1.0) - state.ν[i, j]) - γv * ν_star
end
