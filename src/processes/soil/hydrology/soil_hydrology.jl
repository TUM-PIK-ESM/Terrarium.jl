"""
Base type for implementations of soil water flow dynamics.
"""
abstract type AbstractSoilWaterFluxOperator <: AbstractOperator end

"""
Represents the simplest case of immobile soil water.
"""
struct NoFlow <: AbstractSoilWaterFluxOperator end

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
struct SoilHydrology{
    NF,
    Operator<:AbstractSoilWaterFluxOperator,
    SoilHydraulics<:AbstractSoilHydraulics{NF}
} <: AbstractSoilHydrology{NF}
    "Soil water flux operator"
    operator::Operator

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulics
end

SoilHydrology(
    ::Type{NF},
    operator::AbstractSoilWaterFluxOperator = NoFlow();
    hydraulic_properties::AbstractSoilHydraulics = SoilHydraulicsSURFEX(NF),
) where {NF} = SoilHydrology(operator, hydraulic_properties)

"""
    get_swrc(hydrology::SoilHydrology)

Return the soil water retention curve from the `hydraulic_properties` associated with
the given `SoilHydrology` configuration.
"""
get_swrc(hydrology::SoilHydrology) = hydrology.hydraulic_properties.cond_unsat.swrc

"""
    porosity(idx, state, hydrology::SoilHydrology, strat::AbstractStratigraphy, bgc::AbstractSoilBiogeochemistry)

Return the porosity of the soil volume at `idx` given the current state, hydrology, stratigraphy, and biogeochemistry configurations.
"""
@inline function porosity(idx, state, hydrology::SoilHydrology, strat::AbstractStratigraphy, bgc::AbstractSoilBiogeochemistry)
    org = organic_fraction(idx, state, bgc)
    texture = soil_texture(idx, state, strat)
    return (1 - org)*mineral_porosity(hydrology.hydraulic_properties, texture) + org*organic_porosity(idx, state, bgc)
end

# Immobile soil water (NoFlow)

variables(::SoilHydrology{NF,NoFlow}) where {NF} = (
    auxiliary(:saturation_water_ice, XYZ(), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
)

@inline compute_auxiliary!(state, model, hydrology::SoilHydrology) = nothing

@inline compute_tendencies!(state, model, strat::SoilHydrology{NF, NoFlow}) where {NF} = nothing
