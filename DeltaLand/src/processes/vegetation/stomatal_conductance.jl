"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
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

function compute_auxiliary!(state, stomcond::MedlynStomatalConductance)
    grid = get_grid(photo)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, stomcond)
end

@kernel function compute_auxiliary_kernel!(state, stomcond::MedlynStomatalConductance{NF}) where NF
    i, j = @index(Global, NTuple)

    # Compute Vapor Pressure Deficit [Pa]
    vpd = compute_vpd(state.T_air[i, j], state.q_air[i, j], state.pres[i, j])

    # Compute λc, derived from the optimal stomatal conductance model (Medlyn et al. 2011)
    state.λc[i, j] = NF(1.0) - NF(1.6) / (NF(1.0) + stomcond.g1 / sqrt(vpd))
end
