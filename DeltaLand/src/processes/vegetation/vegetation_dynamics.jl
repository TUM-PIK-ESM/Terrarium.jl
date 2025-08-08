"""
    $TYPEDEF

Vegetation dynamics implementation following PALADYN (Willeit 2016) for a single PFT
based on the Lotka–Volterra approach.

Authors: Maha Badri 

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

"""
    $SIGNATURES

Computes the disturbance rate`γv`.
"""
@inline function compute_γv(veg_dynamics::PALADYNVegetationDynamics)
    # TODO add PALADYN implemetation for the disturbance rate (depends on soil moisture)
    # Placeholder for now γv = min. disturbance rate
    return veg_dynamics.γv_min
end

"""
    $SIGNATURES

Computes `ν_star` which is the maximum between the current vegetation fraction `ν` and the seed fraction `ν_seed`,
to ensure that a PFT is always seeded.
"""
@inline function compute_ν_star(veg_dynamics::PALADYNVegetationDynamics, ν) 
    return max(ν, veg_dynamics.ν_seed)
end


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
    γv = compute_γv(carbon_dynamics)

    # Compute ν_star
    ν_star = compute_ν_star(veg_dynamics, state.ν[i, j])

    # Compute the vegetation fraction tendency
    state.ν_tendency[i, j] = (λ_NPP * state.NPP[i, j] / state.C_veg[i, j]) * 
                             ν_star * (NF(1.0) - state.ν[i, j]) - γv * ν_star
end
