abstract type AbstractSoilHydraulicProperties end

"""
    $TYPEDEF

SURFEX parameterization of mineral soil porosity (Masson et al. 2013).
"""
@kwdef struct SURFEXHydraulics{NF} <: AbstractSoilHydraulicProperties
    porosity::NF = 0.49
    porosity_sand_coef::NF = -1.1e-3
    wilting_point_coef::NF = 37.1e-3
    field_capacity_coef::NF = 89.0e-3
    field_capacity_exp::NF = 0.35
end

# TODO: Maybe we can borrow something better from SINDABD here; the SURFEX scheme is quite simplistic

@inline function mineral_porosity(hydraulics::SURFEXHydraulics, texture::SoilTexture)
    p₀ = hydraulics.porosity
    β = hydraulics.porosity_sand_coef
    por = p₀ + β*texture.sand
    return por
end

@inline function mineral_wilting_point(hydraulics::SURFEXHydraulics, texture::SoilTexture)
    β = hydraulics.wilting_point_coef
    wp = β*sqrt(texture.clay)
    return wp
end

@inline function mineral_field_capacity(hydraulics::SURFEXHydraulics, texture::SoilTexture)
    η = hydraulics.field_capacity_exp
    β = hydraulics.field_capacity_coef
    fc = β*texture.clay^η
    return fc
end
