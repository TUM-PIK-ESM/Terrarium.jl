"""
    $TYPEDEF
Stomatal conductance implementation from PALADYN (Willeit 2016) following the optimal stomatal conductance model
(Medlyn et al. 2011).

Authors: Maha Badri and Matteo Willeit

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
    auxiliary(:T_air, XY()), # Surface air temperature in Celsius [°C]
    auxiliary(:q_air, XY()), # Surface air specific humidity [kg/kg]
    auxiliary(:pres, XY()), # Surface pressure [Pa]
    auxiliary(:λc, XY()), # Ratio of leaf-internal and air CO2 concentration 
)

"""
    $SIGNATURES

Computes the ratio of leaf-internal and air CO2 concentration `λc`, 
derived from the optimal stomatal conductance model (Medlyn et al. 2011).
"""
@inline function compute_λc(stomcond::MedlynStomatalConductance{NF}, vpd) where NF
    λc = NF(1.0) - NF(1.6) / (NF(1.0) + stomcond.g1 / sqrt(vpd))
    return λc
end

function compute_auxiliary!(state, model, stomcond::MedlynStomatalConductance)
    grid = get_grid(model)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, stomcond)
end

@kernel function compute_auxiliary_kernel!(state, stomcond::MedlynStomatalConductance{NF}) where NF
    i, j = @index(Global, NTuple)

    # Compute vpd [Pa]
    vpd = compute_vpd(state.T_air[i, j], state.q_air[i, j], state.pres[i, j])

    # Compute λc
    state.λc[i, j] = compute_λc(stomcond, vpd)
end
