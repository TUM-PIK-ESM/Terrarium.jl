@kwdef struct MedlynStomatalConductance{NF} <: AbstractStomatalConductance
    # TODO check pysical meaning of this parameter
    "Parameter in optimal stomatal conductance formulation (Lin et al. 2015), PFT specific"
    g1::NF = 2.3 # Value for Needleleaf tree PFT 
end

variables(stomcond::MedlynStomatalConductance) = (
    # TODO for now define atmospheric inputs/forcings here, move later
    auxiliary(:T_air, XY()), # Surface air temperature in Kelvin [K]
    auxiliary(:q_air, XY()), # Surface air specific humidity [kg/kg]
    auxiliary(:p, XY()), # Surface pressure [Pa]
    auxiliary(:λc, XY()), # Ratio of leaf-internal and air CO2 concentration 
)

# TODO for now define functions that compute derived variables from atm. inputs/forcings here, move later
@inline function compute_vpd(idx, state, model::AbstractVegetationModel, stomcond::MedlynStomatalConductance{NF}) where NF
    i, j = idx
    
    # Get atmospheric inputs/forcings 
    T_air = state.T_air[i, j] 
    q_air = state.q_air[i, j] 
    p = state.p[i, j] 

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

    # Compute Vapor Pressure Deficit [Pa]
    vpd = compute_vpd(idx, state, model, stomcond)

    # Compute λc, derived from the optimal stomatal conductance model (Medlyn et al. 2011)
    state.λc[i, j] = NF(1.0) - NF(1.6) / (NF(1.0) + stomcond.g1/sqrt(vpd))
end
   