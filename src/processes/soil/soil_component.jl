@kwdef struct Soil{
    NF,
    Stratigraphy <: AbstractStratigraphy,
    Energy <: AbstractSoilEnergyBalance,
    Hydrology <: AbstractSoilHydrology,
    Biogeochemistry <: AbstractSoilBiogeochemistry
} <: AbstractGround{NF}
    "Soil stratigraphy parameterization"
    strat::Stratigraphy

    "Soil energy balance process"
    energy::Energy

    "Soil hydrology (water balance) process"
    hydrology::Hydrology

    "Soil biogeochemistry process"
    biogeochem::Biogeochemistry
end

# Process interface methods

function compute_auxiliary!(
    state, grid,
    soil::Soil,
    constants::PhysicalConstants
)
    # TODO: consider implementing fused kernel here?
    compute_auxiliary!(state, grid, soil.hydrology, soil, constants)
    compute_auxiliary!(state, grid, soil.biogeochem, soil, constants)
    compute_auxiliary!(state, grid, soil.energy, soil, constants)
    return nothing
end

function compute_tendencies!(
    state, grid,
    soil::Soil,
    constants::PhysicalConstants
)
    # TODO: consider implementing fused kernel here?
    compute_tendencies!(state, grid, soil.hydrology, soil, constants)
    compute_tendencies!(state, grid, soil.biogeochem, soil, constants)
    compute_tendencies!(state, grid, soil.energy, soil, constants)
    return nothing
end

# Closures

function closure!(
    state, grid,
    soil::Soil,
    constants::PhysicalConstants
)
    closure!(state, grid, soil.hydrology, soil, constants)
    closure!(state, grid, soil.energy, soil, constants)
end

function invclosure!(
    state, grid,
    soil::Soil,
    constants::PhysicalConstants
)
    invclosure!(state, grid, soil.hydrology, soil, constants)
    invclosure!(state, grid, soil.energy, soil, constants)
end
