@kwdef struct PaladynVegetationDynamics{NF} <: AbstractVegetationDynamics
    "Vegetation seed fraction"
    ν_seed::NF = 0.001

    "Minimum vegetation disturbance rate [1/year]"
    γv_min::NF = 0.002 
end

variables(veg_dynamics::PaladynVegetationDynamics) = (
    prognostic(:veg_fraction, XY()), # PFT fractional area coverage
)

function compute_auxiliary!(idx, state, model, veg_dynamics::PaladynVegetationDynamics)
    i, j = idx
    return nothing
end

function compute_tendencies!(idx, state, model, veg_dynamics::PaladynVegetationDynamics)
    i, j = idx
    
    # Get state variables 
    NPP = state.NPP[i, j]
    C_veg = state.veg_carbon_pool[i, j]
    
    # TODO is this ok or should be recomputed here?
    λ_NPP = state.λ_NPP[i, j]

    # Compute the disturbance rate
    # TODO add Paladyn implemetation
    # Placeholder for now = min dist. rate
    γv = veg_dynamics.γv_min

    # Compute the vegetation fraction tendency
    ν_star = max(state.veg_fraction[i, j], veg_dynamics.ν_seed) 
    state.veg_fraction_tendency[i, j] = (λ_NPP * NPP / C_veg) * ν_star * 
                                        (1.0 - state.veg_fraction[i, j]) - γv * ν_star
end
