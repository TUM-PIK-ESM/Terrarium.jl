"""
    $TYPEDEF

Vegetation dynamics implementation following PALADYN (Willeit 2016) for a single PFT
based on the Lotka–Volterra approach.

Authors: Maha Badri and Matteo Willeit

Properties:
$TYPEDFIELDS
"""
@kwdef struct PALADYNVegetationDynamics{NF} <: AbstractVegetationDynamics
    "Vegetation seed fraction"
    ν_seed::NF = 0.001

    "Minimum vegetation disturbance rate [1/year]"
    γv_min::NF = 0.002 
end

variables(::PALADYNVegetationDynamics) = (
    prognostic(:ν, XY()), # PFT fractional area coverage
    auxiliary(:C_veg, XY()), # Vegetation carbon pool [kgC/m²]
    auxiliary(:LAI_b, XY()), # Balanced Leaf Area Index 
)

function compute_auxiliary!(state, model, veg_dynamics::PALADYNVegetationDynamics)
    # Nothing needed here for now
    return nothing
end

function compute_tendencies!(state, model, veg_dynamics::PALADYNVegetationDynamics)
    grid = get_grid(model)
    launch!(grid, :xy, compute_tendencies_kernel!, state, veg_dynamics, get_carbon_dynamics(model))
end

@kernel function compute_tendencies_kernel!(
    state,
    veg_dynamics::PALADYNVegetationDynamics{NF},
    carbon_dynamics::PALADYNCarbonDynamics{NF}
) where NF
    i, j = @index(Global, NTuple)
    
    # Compute λ_NPP
    λ_NPP = compute_λ_NPP(carbon_dynamics, state.LAI_b[i, j])

    # Compute the disturbance rate
    # TODO add PALADYN implemetation for dist. rate
    # Placeholder for now γv = min. disturbance rate
    γv = veg_dynamics.γv_min

    # Compute the vegetation fraction tendency
    ν_star = max(state.ν[i, j], veg_dynamics.ν_seed) 
    state.ν_tendency[i, j] = (λ_NPP * state.NPP[i, j] / state.C_veg[i, j]) * 
                             ν_star * (NF(1.0) - state.ν[i, j]) - γv * ν_star
end
