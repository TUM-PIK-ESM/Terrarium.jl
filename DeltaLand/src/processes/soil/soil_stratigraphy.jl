@kwdef struct HomogeneousSoil{FC,NF} <: AbstractStratigraphy
    "Material composition of mineral soil componnet"
    texture::SoilTexture{NF,NF,NF} = SoilTexture(:sand)

    "Fraction of pore space that is saturated with water/ice."
    saturation::NF = 1.0

    "Soil freezing characteristic curve"
    freezecurve::FC = FreezeCurves.FreeWater()
end

variables(::HomogeneousSoil) = ()

# TODO: This method interface assumes a single freeze curve for the whole stratigraphy;
# we should ideally relax this assumption for multi-layer stratigraphies, although
# this probably only matters in permafrost environments.
freezecurve(strat::HomogeneousSoil) = strat.freezecurve

# SURFEX parameterization of mineral soil porosity (Masson et al. 2013)
# Maybe we can borrow something better from SINDABD here; this is quite simplistic
mineral_porosity(strat::HomogeneousSoil) = 0.49 - 0.11*strat.texture.sand

pore_water_ice_saturation(strat::HomogeneousSoil) = strat.saturation

@inline pore_water_ice_saturation(i, j, k, state, model, strat::HomogeneousSoil) = pore_water_ice_saturation(strat)

@inline function porosity(i, j, k, state, model, strat::HomogeneousSoil)
    bgc = get_biogeochemistry(model)
    org = organic_fraction(i, j, k, state, model, bgc)
    # note that this assumes the natural porosity of organic material to be constant;
    # probably not a great assumption?
    return (1 - org)*mineral_porosity(strat) + org*organic_porosity(bgc)
end

# do nothing for now
@inline update_state(i, j, k, state, model, strat::HomogeneousSoil) = nothing

@inline compute_tendencies(i, j, k, state, model, strat::HomogeneousSoil) = nothing
