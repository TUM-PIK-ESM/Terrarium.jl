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
    "Parameter in optimal stomatal conductance formulation, Lin et al. 2015 [-], PFT specific"
    g1::NF = 2.3 # Value for Needleleaf tree PFT 
end

MedlynStomatalConductance(::Type{NF}; kwargs...) where {NF} = MedlynStomatalConductance(; kwargs...)

variables(::MedlynStomatalConductance) = (
    auxiliary(:λc, XY()), # Ratio of leaf-internal and air CO2 concentration [-]
)

"""
    $SIGNATURES

Computes the ratio of leaf-internal and air CO2 concentration `λc`, 
derived from the optimal stomatal conductance model (Medlyn et al. 2011),
Eq. 71, PALADYN (Willeit 2016).
"""
@inline function compute_λc(stomcond::MedlynStomatalConductance{NF}, vpd) where NF
    λc = NF(1.0) - NF(1.6) / (NF(1.0) + stomcond.g1 / sqrt(vpd * NF(1.0e-3)))
    return λc
end

function compute_auxiliary!(state, model, stomcond::MedlynStomatalConductance)
    grid = get_grid(model)
    bcs = get_boundary_conditions(model)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, stomcond, bcs.top)
end

@kernel function compute_auxiliary_kernel!(
    state,
    stomcond::MedlynStomatalConductance{NF},
    ::PrescribedAtmosphere
) where NF
    i, j = @index(Global, NTuple)

    # Get inputs
    T_air = state.T_air[i, j] # °C
    q_air = state.q_air[i, j] # kg/kg
    pres = state.pressure[i, j] # Pa

    # Compute vpd [Pa]
    vpd = compute_vpd(T_air, q_air, pres)

    # Compute λc
    λc = compute_λc(stomcond, vpd)

    # Store result
    state.λc[i, j] = λc
end
