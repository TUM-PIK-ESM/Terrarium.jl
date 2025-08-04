abstract type AbstractSoilWaterFluxes end
struct NoFlow <: AbstractSoilWaterFluxes end

@kwdef struct SoilHydrology{
    SoilWaterFluxes<:AbstractSoilWaterFluxes,
    SoilHydraulicProperties<:AbstractSoilHydraulicProperties,
    FC<:FreezeCurves.FreezeCurve
} <: AbstractSoilHydrology
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

variables(::SoilHydrology{NoFlow}) = (
    auxiliary(:pore_water_ice_saturation, XYZ()),
)

@inline function porosity(idx, state, model::AbstractSoilModel, hydrology::SoilHydrology)
    bgc = get_biogeochemistry(model)
    strat = get_stratigraphy(model)
    props = get_hydraulic_properties(hydrology)
    org = organic_fraction(idx, state, model, bgc)
    # note that this currently assumes the natural porosities of both the mineral and organic components to be constant;
    # we probably want to relax this in the future since there are global products of soil texture at least for the upper layers
    texture = get_soil_texture(strat)
    return (1 - org)*mineral_porosity(props, texture) + org*organic_porosity(bgc)
end

@inline compute_auxiliary!(idx, state, model::AbstractSoilModel, hydrology::SoilHydrology) = nothing

@inline compute_tendencies!(idx, state, model, strat::SoilHydrology{NoFlow}) = nothing
