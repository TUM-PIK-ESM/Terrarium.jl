struct SoilEnergyWaterCarbon{
    NF,
    Stratigraphy <: AbstractStratigraphy{NF},
    Energy <: AbstractSoilEnergyBalance{NF},
    Hydrology <: AbstractSoilHydrology{NF},
    Biogeochemistry <: AbstractSoilBiogeochemistry{NF}
} <: AbstractSoil{NF}
    "Soil stratigraphy parameterization"
    strat::Stratigraphy

    "Soil energy balance process"
    energy::Energy

    "Soil hydrology (water balance) process"
    hydrology::Hydrology

    "Soil biogeochemistry process"
    biogeochem::Biogeochemistry
end

function SoilEnergyWaterCarbon(
    ::Type{NF};
    strat = HomogeneousStratigraphy(NF),
    energy = SoilEnergyBalance(NF),
    hydrology = SoilHydrology(NF),
    biogeochem = ConstantSoilCarbonDensity(NF)
) where {NF}
    return SoilEnergyWaterCarbon(strat, energy, hydrology, biogeochem)
end

# Process interface methods

function initialize!(
    state, grid,
    soil::SoilEnergyWaterCarbon,
    constants::PhysicalConstants
)
    initialize!(state, grid, soil.hydrology, soil, constants)
    initialize!(state, grid, soil.biogeochem, soil, constants)
    initialize!(state, grid, soil.energy, soil, constants)
    return nothing
end

function compute_auxiliary!(
    state, grid,
    soil::SoilEnergyWaterCarbon,
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
    soil::SoilEnergyWaterCarbon,
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
    soil::SoilEnergyWaterCarbon,
    constants::PhysicalConstants
)
    closure!(state, grid, get_closure(soil.hydrology), soil.hydrology, soil)
    closure!(state, grid, get_closure(soil.energy), soil.energy, soil, constants)
end

function invclosure!(
    state, grid,
    soil::SoilEnergyWaterCarbon,
    constants::PhysicalConstants
)
    invclosure!(state, grid, get_closure(soil.hydrology), soil.hydrology, soil)
    invclosure!(state, grid, get_closure(soil.energy), soil.energy, soil, constants)
end
