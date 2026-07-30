"""
    $TYPEDEF

Coupled process type that encapsulates the coupling of soil energy, water, and carbon dynamics.
The stratigraphy parameterization determines how the vertical layering of the soil is parameterized.
"""
struct SoilEnergyWaterCarbon{
        NF,
        Stratigraphy <: AbstractStratigraphy,
        Energy <: AbstractSoilThermodynamics,
        Hydrology <: AbstractSoilHydrology,
        Biogeochemistry <: AbstractSoilBiogeochemistry,
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

# No component is pinned to `NF` — any of them may hold a promoted (e.g. traced) parameter — so `NF`
# is not derivable from the field types and must be given. `constructorof` rebuilds through this
# constructor so `setproperties` / `Flatten.reconstruct` / `Adapt` retain `NF`.
SoilEnergyWaterCarbon{NF}(
    strat::AbstractStratigraphy,
    energy::AbstractSoilThermodynamics,
    hydrology::AbstractSoilHydrology,
    biogeochem::AbstractSoilBiogeochemistry,
) where {NF} = SoilEnergyWaterCarbon{NF, typeof(strat), typeof(energy), typeof(hydrology), typeof(biogeochem)}(
    strat, energy, hydrology, biogeochem
)

ConstructionBase.constructorof(::Type{<:SoilEnergyWaterCarbon{NF}}) where {NF} = SoilEnergyWaterCarbon{NF}

function SoilEnergyWaterCarbon(
        ::Type{NF};
        strat = HomogeneousSoilStratigraphy(NF),
        energy = SoilThermodynamics(NF),
        hydrology = SoilHydrology(NF),
        biogeochem = ConstantSoilCarbonDensity(NF)
    ) where {NF}
    return SoilEnergyWaterCarbon{NF}(strat, energy, hydrology, biogeochem)
end

# Process interface methods

"""
    $TYPEDEF

Initialize the soil energy, water, and carbon state variables on `grid` given
the parameter values in `constants`.
"""
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

"""
    $TYPEDEF

Compute auxiliary variables for soil energy, water, and carbon state variables
on `grid` based on the given values in `constants`.
"""
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

"""
    $TYPEDEF

Compute tendencies for soil energy, water, and carbon state variables on `grid`
based on the given values in `constants`.
"""
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

"""
    $TYPEDEF

Compute the forward closure mapping for soil hydrology and energy, in that order.
"""
function closure!(
        state, grid,
        soil::SoilEnergyWaterCarbon,
        constants::PhysicalConstants
    )
    closure!(state, grid, get_closure(soil.hydrology), soil.hydrology, soil)
    closure!(state, grid, get_closure(soil.energy), soil.energy, soil, constants)
    return nothing
end

"""
    $TYPEDEF

Compute the inverse closure mapping for soil hydrology and energy, in that order.
"""
function invclosure!(
        state, grid,
        soil::SoilEnergyWaterCarbon,
        constants::PhysicalConstants
    )
    invclosure!(state, grid, get_closure(soil.hydrology), soil.hydrology, soil)
    invclosure!(state, grid, get_closure(soil.energy), soil.energy, soil, constants)
    return nothing
end
