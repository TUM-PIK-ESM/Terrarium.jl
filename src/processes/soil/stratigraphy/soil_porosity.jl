"""
    $TYPEDEF

Parameterization of soil porosity that simply specifies constant values
for the `mineral` and `organic` components.
"""
@kwdef struct ConstantSoilPorosity{NF} <: AbstractSoilPorosity{NF}
    "Prescribed mineral soil porosity [-]"
    mineral_porosity::NF = 0.49

    "Natural porosity of organic material"
    organic_porosity::NF = 0.90
end

@inline organic_porosity(props::ConstantSoilPorosity, texture::SoilTexture) = props.organic_porosity

@inline mineral_porosity(props::ConstantSoilPorosity, texture::SoilTexture) = props.mineral_porosity

"""
    $TYPEDEF

SURFEX parameterization of mineral soil porosity (Masson et al. 2013).
"""
@kwdef struct SoilPorositySURFEX{NF} <: AbstractSoilPorosity{NF}
    "Assumed base porosity of the soil without sand [-]"
    porosity_base::NF = 0.49

    "Linear coeficient of porosity adjustment for sand [-]"
    porosity_sand_coef::NF = -0.11

    "Natural porosity of organic material"
    porosity_organic::NF = 0.90
end

@inline organic_porosity(props::SoilPorositySURFEX) = por.porosity_organic

@inline function mineral_porosity(props::SoilPorositySURFEX)
    p₀ = props.porosity_base
    β_s = props.porosity_sand_coef
    por = p₀ + β_s*texture.sand
    return por
end

# Kernel functions

mineral_porosity(i, j, k, grid, state, props::AbstractSoilPorosity, texture::SoilTexture) = mineral_porosity(props, texture)

organic_porosity(i, j, k, grid, state, props::AbstractSoilPorosity, texture::SoilTexture) = organic_porosity(props, texture)
