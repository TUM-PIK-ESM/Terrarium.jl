abstract type AbstractSoilHydraulicProperties end

"""
    $TYPEDEF

SURFEX parameterization of mineral soil porosity (Masson et al. 2013).
"""
@kwdef struct SURFEXSoilHydraulics{NF} <: AbstractSoilHydraulicProperties
    porosity::NF = 0.49
    porosity_sand_coef::NF = -0.11
    wilting_point::NF = 0.05
end

# Maybe we can borrow something better from SINDABD here; the SURFEX scheme is quite simplistic
mineral_porosity(hydraulics::SURFEXSoilHydraulics, texture::SoilTexture) = hydraulics.porosity + hydraulics.*texture.sand

wilting_point(hydraulics::SURFEXSoilHydraulics) = hydraulics.porosity_sand_coef
