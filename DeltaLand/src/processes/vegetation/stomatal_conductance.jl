@kwdef struct MedlynStomatalConductance{NF} <: AbstractStomatalConductance
    # TODO check g1 parameter physical meaning
    "Parameter in optimal stomatal conductance formulation (Lin et al. 2015), PFT specific"
    g1::NF = 2.3 # Value for Needleleaf tree PFT (always evergreen)
end

variables(stomcond::MedlynStomatalConductance) = (
    auxiliary(:λc, XY()), # Ratio of leaf-internal and air CO2 concentration (λc)
)

@inline function compute_vpd(idx, state, model::AbstractVegetationModel, stomcond::MedlynStomatalConductance{NF}) where NF
    # TODO needs T_air_C, p, q_air as input, also not really specific to stomatal conductance?
    # Compute Saturation vapor pressure over water [Pa]
    e_sat_w = NF(6.1094e2) * exp(NF(17.625) * T_air_C / (NF(243.04) + T_air_C))

    # Convert air specific humidity to vapor pressure [Pa]
    q_to_e = q_air * p / (NF(0.622) + NF(0.378) * q_air)

    # Compute vapor pressure deficit [Pa]
    # TODO is the max operation needed?
    vpd = max(e_sat_w - q_to_e, NF(0.1))
    
    return vpd
end

@inline function compute_auxiliary!(idx, state, model::AbstractVegetationModel, stomcond::MedlynStomatalConductance{NF}) where NF
    i, j = idx

    # Compute Vapor Pressure Deficit (Pa)
    vpd = compute_vpd(idx, state, model, stomcond)

    # Compute λc, derived from the optimal stomatal conductance model (Medlyn et al. 2011)
    state.λc[i, j] = NF(1.0) - NF(1.6) / (NF(1.0) + stomcond.g1/sqrt(vpd))
end
   