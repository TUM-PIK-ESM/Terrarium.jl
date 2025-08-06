@kwdef struct MedlynStomatalConductance{NF} <: AbstractStomatalConductance
    # TODO check g1 parameter physical meaning
    "Parameter in optimal stomatal conductance formulation (Lin et al. 2015), PFT specific"
    g1::NF = 2.3 # Value for Needleleaf tree PFT (always evergreen)
end

variables(stomcond::MedlynStomatalConductance) = (
    auxiliary(:λc, XY()), # Ratio of leaf-internal and air CO2 concentration (λc)
)

@inline function compute_auxiliary!(idx, state, model::AbstractVegetationModel, stomcond::MedlynStomatalConductance{NF}) where NF
    i, j = idx

    # Compute Vapor Pressure Deficit (Pa)
    # TODO add VPD implementation
    # For now, set vpd to a random value
    vpd = rand(NF)

    # Compute λc, derived from the optimal stomatal conductance model (Medlyn et al. 2011)
    state.λc[i, j] = NF(1.0) - NF(1.6) / (NF(1.0) + stomcond.g1/sqrt(vpd))
end
   