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
    hydraulic_properties::SoilHydraulicProperties = SURFEXSoilHydraulics()
    
    "Soil freezing characteristic curve"
    freezecurve::FC = FreezeCurves.FreeWater()
end

variables(::SoilHydrology{NoFlow}) = (
    auxiliary(:pore_water_ice_saturation, XYZ())
)

# TODO: This method interface assumes a single freeze curve for the whole stratigraphy;
# we should ideally relax this assumption for multi-layer stratigraphies, although
# this probably only matters in permafrost environments.
freezecurve(hydrology::SoilHydrology) = hydrology.freezecurve

pore_water_ice_saturation(hydrology::SoilHydrology) = hydrology.saturation

@inline pore_water_ice_saturation(idx, state, model, hydrology::SoilHydrology) = pore_water_ice_saturation(hydrology)

@inline function porosity(idx, state, model, hydrology::SoilHydrology)
    bgc = get_biogeochemistry(model)
    strat = get_stratigraphy(model)
    org = organic_fraction(idx, state, model, bgc)
    # note that this currently assumes the natural porosities of both the mineral and organic components to be constant;
    # we probably want to relax this in the future since there are global products of soil texture at least for the upper layers
    texture = get_soil_texture(strat)
    return (1 - org)*mineral_porosity(hydrology, texture) + org*organic_porosity(bgc)
end

@inline function update_state(idx, state, model::AbstractSoilModel, hydrology::SoilHydrology)
    i, j, k = idx
    hydrology = get_hydrology(model)
    # set saturation level of pore water/ice to value specified by stratigraphy
    state.pore_water_ice_saturation[i, j, k] = pore_water_ice_saturation(idx, state, model, hydrology)
end

@inline compute_tendencies(idx, state, model, strat::SoilHydrology{NoFlow}) = nothing
