"""
    $TYPEDEF

Vegetation dynamics implementation following PALADYN (Willeit 2016) for a single PFT
based on the Lotka–Volterra approach.

Authors: Maha Badri 

Properties:
$TYPEDFIELDS
"""
@kwdef struct PALADYNVegetationDynamics{NF} <: AbstractVegetationDynamics
    "Vegetation seed fraction [-]"
    ν_seed::NF = 0.001

    "Minimum vegetation disturbance rate [1/year]"
    # TODO this parameter is yearly, should be changed to daily for now
    γv_min::NF = 0.002 
end

variables(::PALADYNVegetationDynamics) = (
    prognostic(:ν, XY()), # PFT fractional area coverage [-]
    auxiliary(:C_veg, XY()), # Vegetation carbon pool [kgC/m²]
    auxiliary(:LAI_b, XY()), # Balanced Leaf Area Index [m²/m²]
)

"""
    $SIGNATURES

Computes the disturbance rate`γv`,
Eq. 80, PALADYN (Willeit 2016).
"""
@inline function compute_γv(veg_dynamics::PALADYNVegetationDynamics)
    # TODO add PALADYN implemetation for the disturbance rate (depends on soil moisture)
    # Placeholder for now γv = min. disturbance rate
    return veg_dynamics.γv_min
end

"""
    $SIGNATURES

Computes `ν_star` which is the maximum between the current vegetation fraction `ν` and the seed fraction `ν_seed` [-],
to ensure that a PFT is always seeded.
"""
@inline function compute_ν_star(veg_dynamics::PALADYNVegetationDynamics, ν) 
    return max(ν, veg_dynamics.ν_seed)
end

"""
    $SIGNATURES

Computes the vegetation fraction tendency for a single PFT,
Eq. 73, PALADYN (Willeit 2016).
"""
@inline function compute_ν_tend(
    veg_dynamics::PALADYNVegetationDynamics, 
    vegcarbon_dynamics::PALADYNCarbonDynamics{NF},
    LAI_b::NF,
    C_veg::NF,
    ν::NF) where NF

    # Compute λ_NPP
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b)

    # Compute the disturbance rate
    γv = compute_γv(veg_dynamics)

    # Compute ν_star
    ν_star = compute_ν_star(veg_dynamics, ν)

    # Compute the vegetation fraction tendency 
    ν_tendency =  (λ_NPP * C_veg / ν_star) * (NF(1.0) - ν) - γv * ν_star
    return ν_tendency
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
    vegcarbon_dynamics::PALADYNCarbonDynamics{NF}
) where NF
    i, j = @index(Global, NTuple)

    # Get inputs
    LAI_b = state.LAI_b[i, j]
    C_veg = state.C_veg[i, j]
    ν = state.ν[i, j]

    # Compute the vegetation fraction tendency
    ν_tendency = compute_ν_tend(veg_dynamics, vegcarbon_dynamics, LAI_b, C_veg, ν)
    
    # Store result
    state.ν_tendency[i, j] = ν_tendency
end
