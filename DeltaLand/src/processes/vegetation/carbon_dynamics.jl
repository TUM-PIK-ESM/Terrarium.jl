@kwdef struct PaladynCarbonDynamics{NF} <: AbstractVegetationCarbonDynamics
    "Specific leaf area (Kattge et al. 2011) [m²/kgC], PFT specific"
    SLA::NF = 10 # Value for Needleleaf tree PFT (always evergreen)
    
    "Allometric coefficient [kgC/m²], PFT specific"
    # TODO check awl value in Paladyn code
    awl::NF = 2.0 # Value for Needleleaf tree PFT (always evergreen) 
    
    "Minimum Leaf Area Index modified from Clark et al. 2011, PFT specific"
    LAI_min::NF = 1.0 # Value for Needleleaf tree PFT (always evergreen)
    
    "Maximum Leaf Area Index modified from Clark et al. 2011, PFT specific"
    LAI_max::NF = 6.0 # Value for Needleleaf tree PFT (always evergreen)
    
    "Leaf turnover rate (Kattge et al. 2011) [1/year], PFT specific"
    γL ::NF = 0.3 # Value for Needleleaf tree PFT (always evergreen)
    
    "Root turnover rate [1/year], PFT specific"
    γR::NF = 0.3 # Value for Needleleaf tree PFT (always evergreen)
    
    "Stem turnover rate modified from Clark et al. 2011 [1/year], PFT specific"
    γS::NF = 0.05 # Value for Needleleaf tree PFT (always evergreen)
end

variables(vegcarbon_dynamics::PaladynCarbonDynamics) = (
    prognostic(:C_veg, XY()), # Vegetation carbon pool (kgC/m²)
    auxiliary(:LAI_b, XY()), # Balanced Leaf Area Index (LAI_b) 
    auxiliary(:NPP, XY()), # Net Primary Production (NPP) kgC/m²/day
)

@inline function compute_auxiliary!(idx, state, model::AbstractVegetationModel, vegcarbon_dynamics::PaladynCarbonDynamics)
    i, j = idx

    # Compute balanced Leaf Area Index (LAI_b)
    # Following Paladyn approach (assuming with bwl = 1)
    state.LAI_b[i, j] = ((2.0 / vegcarbon_dynamics.SLA) + vegcarbon_dynamics.awl) / state.veg_carbon_pool[i, j]

    # TODO move to a separate function
    # Compute λ_NPP based on the balanced Leaf Area Index (LAI_b)
    LAI_b = state.LAI_b[i, j]
    if LAI_b < vegcarbon_dynamics.LAI_min
        state.λ_NPP[i, j] = 0.0
    elseif LAI_b <= vegcarbon_dynamics.LAI_max
        state.λ_NPP[i, j] = (LAI_b - vegcarbon_dynamics.LAI_min) /  
                (vegcarbon_dynamics.LAI_max - vegcarbon_dynamics.LAI_min)
    else
        state.λ_NPP[i, j] = 1.0
    end
end

@inline function compute_tendencies!(idx, state, model::AbstractVegetationModel, vegcarbon_dynamics::PaladynCarbonDynamics)
    i, j = idx

    # Get state variables
    LAI_b = state.LAI_b[i, j]

    # TODO call a separate function to compute λ_NPP
    λ_NPP = state.λ_NPP[i, j]

    # Compute local litterfall rate
    Λ_loc = (vegcarbon_dynamics.γL/vegcarbon_dynamics.SLA +
             vegcarbon_dynamics.γR/vegcarbon_dynamics.SLA + 
             vegcarbon_dynamics.γS*vegcarbon_dynamics.awl) * LAI_b
    
    # Compute the vegetation carbon pool tendency
    state.C_veg_tendency[i, j] = (1.0 - λ_NPP) * state.NPP[i, j] - Λ_loc
end
