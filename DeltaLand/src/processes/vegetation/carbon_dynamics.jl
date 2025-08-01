struct VegetationCarbonDynamics <: AbstractVegetationCarbonDynamics

end

variables(carbon_dynamics::VegetationCarbonDynamics) = (
    prognostic(:carbon_pool, XY()),

)

function compute_tendencies!(idx, state, model, carbon_dynamics::VegetationCarbonDynamics)
    GPP = state.GPP
end
