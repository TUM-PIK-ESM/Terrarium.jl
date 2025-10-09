import FreezeCurves: BrooksCorey, VanGenuchten

"""
Base type for implementations of soil water flow dynamics.
"""
abstract type AbstractSoilWaterFlowOperator <: AbstractOperator end

"""
Represents the simplest case of immobile soil water.
"""
struct NoFlow <: AbstractSoilWaterFlowOperator end

"""
    RichardsEq{NF} <: AbstractSoilWaterFlowOperator

Operator for soil hydrology corresponding to the Richardson-Richards equation for variably saturated
flow in porous media.
"""
@kwdef struct RichardsEq{NF} <: AbstractSoilWaterFlowOperator
    "Exponential scaling factor for ice impedance"
    Ω::NF = 7

    "Closure relation for mapping between water potential (hydraulic head) and saturation"
    saturation_closure = PressureSaturationRelation()
end

get_closure(op::RichardsEq) = op.saturation_closure

struct PressureSaturationRelation <: AbstractClosureRelation end

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
struct SoilHydrology{
    NF,
    Operator<:AbstractSoilWaterFlowOperator,
    RetentionCurve<:SWRC,
    SoilHydraulicProperties<:AbstractSoilHydraulicProperties{NF}
} <: AbstractSoilHydrology{NF}
    "Soil water flow scheme"
    operator::Operator

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulicProperties
    
    "Soil water retention curve(s) from FreezeCurves.jl"
    swrc::RetentionCurve
end

SoilHydrology(
    ::Type{NF};
    operator::AbstractSoilWaterFlowOperator = NoFlow(),
    hydraulic_properties::AbstractSoilHydraulicProperties{NF} = SURFEXHydraulics(NF),
    freezecurve::FreezeCurve = default_swrc(flow, hydraulic_properties)
) where {NF} = SoilHydrology(operator, hydraulic_properties, freezecurve)

"""
    default_swrc(::AbstractSoilWaterFlowOperator, ::AbstractSoilHydraulicProperties)

Return the default soil water retention curve (SWRC) for the given soil hydrology configuration.
Defaults to `nothing` which represents no use of a pressure-saturation relation.
"""
default_swrc(::AbstractSoilWaterFlowOperator, ::AbstractSoilHydraulicProperties) = nothing

get_hydraulic_properties(hydrology::SoilHydrology) = hydrology.hydraulic_properties

# TODO: This method interface assumes a single water retenction curve for the whole stratigraphy;
# we should ideally relax this assumption for multi-layer stratigraphies
get_soil_water_retention_curve(hydrology::SoilHydrology) = hydrology.swrc

"""
    porosity(idx, state, hydrology::SoilHydrology, strat::AbstractStratigraphy, bgc::AbstractSoilBiogeochemistry)

Return the porosity of the soil volume at `idx` given the current state, hydrology, stratigraphy, and biogeochemistry configurations.
"""
@inline function porosity(idx, state, hydrology::SoilHydrology, strat::AbstractStratigraphy, bgc::AbstractSoilBiogeochemistry)
    props = get_hydraulic_properties(hydrology)
    org = organic_fraction(idx, state, bgc)
    texture = soil_texture(idx, state, strat)
    return (1 - org)*mineral_porosity(props, texture) + org*organic_porosity(idx, state, bgc)
end

# Immobile soil water (NoFlow)

variables(::SoilHydrology{NF,NoFlow}) where {NF} = (
    auxiliary(:saturation_water_ice, XYZ(), bounds=0..1, desc="Saturation level of water and ice in the pore space"),
)

@inline compute_auxiliary!(state, model, hydrology::SoilHydrology) = nothing

@inline compute_tendencies!(state, model, strat::SoilHydrology{NF,NoFlow}) where {NF} = nothing

# TODO: Richardson-Richards equation diffusion/advection

variables(hydrology::SoilHydrology{NF,<:RichardsEq}) where {NF} = (
    prognosic(:matric_potential, XYZ()),
    auxiliary(:saturation_water_ice, XYZ(), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
)
