"""
Base type for implementations of soil water flow dynamics.
"""
abstract type AbstractVerticalFlow <: AbstractOperator end

"""
Represents a hydrology scheme where soil water is immobile.
"""
struct NoFlow <: AbstractVerticalFlow end

"""
Base type for soil/subsurface runoff schemes.
"""
abstract type AbstractSoilRunoff end

"""
Base type for evapotranspirative forcing terms in soil layers.
"""
abstract type AbstractSoilET end

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
struct SoilHydrology{
    NF,
    VerticalFlow<:AbstractVerticalFlow,
    Runoff<:AbstractSoilRunoff,
    SoilET<:Union{Nothing, AbstractSoilET},
    SoilHydraulics<:AbstractSoilHydraulics{NF}
} <: AbstractSoilHydrology{NF}
    "Soil water vertical flow operator"
    vertflow::VerticalFlow

    "Soil subsurface runoff scheme"
    runoff::Runoff

    "Soil evapotranspiration scheme"
    evapotranspiration::SoilET

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulics
end

SoilHydrology(
    ::Type{NF},
    vertflow::AbstractVerticalFlow = NoFlow();
    runoff::AbstractSoilRunoff = SoilRunoff(),
    evapotranspiration::Union{Nothing, AbstractSoilET} = nothing,
    hydraulic_properties::AbstractSoilHydraulics = SoilHydraulicsSURFEX(NF),
) where {NF} = SoilHydrology(vertflow, runoff, evapotranspiration, hydraulic_properties)

"""
    get_swrc(hydrology::SoilHydrology)

Return the soil water retention curve from the `hydraulic_properties` associated with
the given `SoilHydrology` configuration.
"""
get_swrc(hydrology::SoilHydrology) = hydrology.hydraulic_properties.cond_unsat.swrc

"""
State variables for `SoilHydrology` processes.
"""
variables(hydrology::SoilHydrology{NF}) where {NF} = (
    variables(hydrology.vertflow)...,
    variables(hydrology.runoff)...,
    variables(hydrology.evapotranspiration)...,
)

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

variables(::NoFlow) = (
    auxiliary(:saturation_water_ice, XYZ(), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
)

@inline initialize!(state, model, hydrology::SoilHydrology) = nothing

@inline compute_auxiliary!(state, model, hydrology::SoilHydrology{NF, NoFlow}) where {NF} = nothing

@inline compute_tendencies!(state, model, hydrology::SoilHydrology{NF, NoFlow}) where {NF} = nothing

# Runoff

"""
    $TYPEDEF

Generic scheme for respresenting soil/subsurface runoff as a composition of three components:
excess infiltration, interflow, and drainage.

Not yet implemented.
"""
@kwdef struct SoilRunoff{SR, IF, DR} <: AbstractSoilRunoff
    "Excess infiltration runoff scheme"
    excess_infiltration::SR = nothing
    
    "Subsurface interflow runoff"
    interflow::IF = nothing

    "Groundwater drainage runoff"
    drainage::DR = nothing
end
