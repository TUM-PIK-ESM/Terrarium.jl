"""
    $TYPEDEF

Vegetation carbon dynamics implementation following PALADYN (Willeit 2016) but considering only the sum of the vegetation
carbon pools. The subsequent splitting into C_leaf, C_stem, C_root is not implemented for now.

Authors: Maha Badri 

Properties:
$TYPEDFIELDS
"""
@kwdef struct PALADYNCarbonDynamics{NF} <: AbstractVegetationCarbonDynamics
    "Specific leaf area (Kattge et al. 2011) [m²/kgC], PFT specific"
    SLA::NF = 10.0 # Value for Needleleaf tree PFT 
    
    "Allometric coefficient, modified from Cox 2001 to account for bwl=1 [kgC/m²], PFT specific"
    awl::NF = 2.0 # Value for Needleleaf tree PFT 
    
    "Minimum Leaf Area Index modified from Clark et al. 2011, PFT specific"
    LAI_min::NF = 1.0 # Value for Needleleaf tree PFT 
    
    "Maximum Leaf Area Index modified from Clark et al. 2011, PFT specific"
    LAI_max::NF = 6.0 # Value for Needleleaf tree PFT 
    
    "Leaf turnover rate (Kattge et al. 2011) [1/year], PFT specific"
    γL::NF = 0.3 # Value for Needleleaf tree PFT 
    
    "Root turnover rate [1/year], PFT specific"
    γR::NF = 0.3 # Value for Needleleaf tree PFT 
    
    "Stem turnover rate modified from Clark et al. 2011 [1/year], PFT specific"
    γS::NF = 0.05 # Value for Needleleaf tree PFT 
end

variables(::PALADYNCarbonDynamics) = (
    prognostic(:C_veg, XY()), # Vegetation carbon pool [kgC/m²]
    auxiliary(:LAI_b, XY()), # Balanced Leaf Area Index [m²/m²]
    auxiliary(:NPP, XY()), # Net Primary Production [kgC/m²/day]
)

"""
    $SIGNATURES

Computes `λ_NPP`,a factor determining the partitioning of NPP between increase of 
vegetation carbon of the existing vegetated area and spreading of the given PFT
based on the balanced Leaf Area Index `LAI_b`.
"""
@inline function compute_λ_NPP(vegcarbon_dynamics::PALADYNCarbonDynamics{NF}, LAI_b) where NF
    if LAI_b < vegcarbon_dynamics.LAI_min
        λ_NPP = zero(NF)
    elseif LAI_b <= vegcarbon_dynamics.LAI_max
        λ_NPP = (LAI_b - vegcarbon_dynamics.LAI_min) /  
                (vegcarbon_dynamics.LAI_max - vegcarbon_dynamics.LAI_min)
    else
        λ_NPP = NF(1.0)
    end
    return λ_NPP
end

"""
    $SIGNATURES

Computes `LAI_b`, the balanced Leaf Area Index based on the vegetation carbon pool `C_veg`
(assuming with bwl = 1).
"""

@inline function compute_LAI_b(vegcarbon_dynamics::PALADYNCarbonDynamics{NF}, C_veg) where NF   
    LAI_b = ((NF(2.0) / vegcarbon_dynamics.SLA) + vegcarbon_dynamics.awl) / C_veg
    return LAI_b
end

"""
    $SIGNATURES
Computes the local litterfall rate `Λ_loc` based on the balanced Leaf Area Index `LAI_b`
(assuming evergreen PFTs).
"""
@inline function compute_Λ_loc(vegcarbon_dynamics::PALADYNCarbonDynamics{NF}, LAI_b) where NF
    Λ_loc = (vegcarbon_dynamics.γL / vegcarbon_dynamics.SLA +
             vegcarbon_dynamics.γR / vegcarbon_dynamics.SLA + 
             vegcarbon_dynamics.γS * vegcarbon_dynamics.awl) * LAI_b
    return Λ_loc
end
    
function compute_auxiliary!(state, model, vegcarbon_dynamics::PALADYNCarbonDynamics)
    grid = get_grid(model)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, vegcarbon_dynamics)
end

@kernel function compute_auxiliary_kernel!(state, vegcarbon_dynamics::PALADYNCarbonDynamics{NF}) where NF
    i, j = @index(Global, NTuple)

    # Compute balanced Leaf Area Index 
    state.LAI_b[i, j] = compute_LAI_b(vegcarbon_dynamics, state.C_veg[i, j])

end

function compute_tendencies!(state, model, vegcarbon_dynamics::PALADYNCarbonDynamics)
    grid = get_grid(model)
    launch!(grid, :xy, compute_tendencies_kernel!, state, vegcarbon_dynamics)
end

@kernel function compute_tendencies_kernel!(state, vegcarbon_dynamics::PALADYNCarbonDynamics{NF}) where NF  
    i, j = @index(Global, NTuple)

    # Compute λ_NPP
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, state.LAI_b[i, j])

    # Compute local litterfall rate
    Λ_loc = compute_Λ_loc(vegcarbon_dynamics, state.LAI_b[i, j])
    
    # Compute the vegetation carbon pool tendency
    state.C_veg_tendency[i, j] = (NF(1.0) - λ_NPP) * state.NPP[i, j] - Λ_loc
end
