abstract type AbstractSoilWaterFluxes end
struct NoFlow <: AbstractSoilWaterFluxes end
struct RichardsEq{Advection<:AbstractAdvectionScheme} <: AbstractSoilWaterFluxes
    advection::Advection
end

@kwdef struct SoilHydrology{
    NF,
    SoilWaterFluxes<:AbstractSoilWaterFluxes,
    SoilHydraulicProperties<:AbstractSoilHydraulicProperties{NF},
    FC<:FreezeCurves.FreezeCurve
} <: AbstractSoilHydrology{NF}
    "Soil water flux scheme"
    fluxes::SoilWaterFluxes = NoFlow()

    "Soil hydraulic properties parameterization"
    hydraulic_properties::SoilHydraulicProperties = SURFEXHydraulics()
    
    "Soil freezing and water retention curve(s) from FreezeCurves.jl"
    freezecurve::FC = FreezeCurves.FreeWater()
end

get_hydraulic_properties(hydrology::SoilHydrology) = hydrology.hydraulic_properties

# TODO: This method interface assumes a single freeze curve for the whole stratigraphy;
# we should ideally relax this assumption for multi-layer stratigraphies, although
# this probably only matters in permafrost environments.
get_freezecurve(hydrology::SoilHydrology) = hydrology.freezecurve

@inline function porosity(idx, state, hydrology::SoilHydrology, strat::AbstractStratigraphy, bgc::AbstractSoilBiogeochemistry)
    props = get_hydraulic_properties(hydrology)
    org = organic_fraction(idx, state, bgc)
    # note that this currently assumes the natural porosities of both the mineral and organic components to be constant;
    # we probably want to relax this in the future since there are global products of soil texture at least for the upper layers
    texture = get_soil_texture(strat)
    return (1 - org)*mineral_porosity(props, texture) + org*organic_porosity(bgc)
end

# Immobile soil water (NoFlow)

variables(::SoilHydrology{NF,NoFlow}) where {NF} = (
    auxiliary(:pore_water_ice_saturation, XYZ()),
)

@inline compute_auxiliary!(state, model, hydrology::SoilHydrology) = nothing

@inline compute_tendencies!(state, model, strat::SoilHydrology{NF,NoFlow}) where {NF} = nothing

# Richardson-Richards equation diffusion/advection

variables(::SoilHydrology{NF,<:RichardsEq}) where {NF} = (
    prognosic(:pore_water_ice_potential, XYZ()),
    auxiliary(:pore_water_ice_saturation, XYZ()),
)