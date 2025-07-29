@kwdef struct HomogeneousSoil{FC,NF} <: AbstractStratigraphy
    "Material composition of mineral soil componnet"
    texture::SoilTexture{NF,NF,NF} = SoilTexture(:sand)

    "Soil freezing characteristic curve"
    freezecurve::FC = FreeWater()
end

# SURFEX parameterization of mineral soil porosity (Masson et al. 2013)
# Maybe we can borrow something better from SINDABD here; this is quite simplistic
mineral_porosity(soil::HomogeneousSoil) = 0.49 - 0.11*soil.texture.sand

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
