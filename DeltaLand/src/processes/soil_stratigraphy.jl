@kwdef struct HomogeneousSoil{FC,NF} <: AbstractStratigraphy
    "Material composition of mineral soil componnet"
    texture::SoilTexture{NF,NF,NF} = SoilTexture(:sand)

    "Natural porosity of organic material"
    por_org::NF = 0.90 # note: maybe this should be handled by the biogeochem model?

    "Soil freezing characteristic curve"
    freezecurve::FC = FreeWater()
end

# do nothing for now
@inline update_state(i, j, k, state, model, strat::HomogeneousSoil) = nothing

@inline compute_tendencies(i, j, k, state, model, strat::HomogeneousSoil) = nothing
