abstract type AbstractSoilHydraulicProperties{NF} end

"""
    PrescribedHydraulics{NF} <: AbstractSoilHydraulicProperties

Represents a simple case where soil hydraulic properties are given as constant values.
This is mostly provided just for testing, although it may be useful in certain cases where direct
measurements of hydraulic properites are available.

Properties:
$TYPEDFIELDS
"""
@kwdef struct PrescribedHydraulics{NF} <: AbstractSoilHydraulicProperties{NF}
    "Hydraulic conductivity at saturation [m/s]"
    cond_sat::NF = 1e-5

    "Prescribed soil porosity [-]"
    porosity::NF = 0.49

    "Prescribed field capacity [-]"
    field_capacity::NF = 0.1

    "Prescribed wilting point [-]"
    wilting_point::NF = 0.05
end

PrescribedHydraulics(::Type{NF}; kwargs...) where {NF} = PrescribedHydraulics{NF}(; kwargs...)

@inline saturated_hydraulic_conductivity(hydraulics::PrescribedHydraulics, args...) = hydraulics.cond_sat

@inline mineral_porosity(hydraulics::PrescribedHydraulics, args...) = hydraulics.porosity

@inline mineral_wilting_point(hydraulics::PrescribedHydraulics, args...) = hydraulics.wilting_point

@inline mineral_field_capacity(hydraulics::PrescribedHydraulics, args...) = hydraulics.field_capacity

"""
    $TYPEDEF

SURFEX parameterization of mineral soil porosity (Masson et al. 2013).
"""
@kwdef struct SURFEXHydraulics{NF} <: AbstractSoilHydraulicProperties{NF}
    "Hydraulic conductivity at saturation [m/s]"
    cond_sat::NF = 1e-5

    "Base porosity of soil without any sand [-]"
    porosity::NF = 0.49

    "Linear coeficient of porosity adjustment due to sand content [-]"
    porosity_sand_coef::NF = -1.1e-3
    
    "Linear coeficient of wilting point adjustment due to clay content [-]"
    wilting_point_coef::NF = 37.13e-3
    
    "Linear coeficient of field capacity adjustment due to clay content [-]"
    field_capacity_coef::NF = 89.0e-3
    
    "Exponent of field capacity adjustment due to clay content [-]"
    field_capacity_exp::NF = 0.35
end

SURFEXHydraulics(::Type{NF}; kwargs...) where {NF} = SURFEXHydraulics{NF}(; kwargs...)

# TODO: this is not quite correct, SURFEX uses a hydraulic conductivity function that decreases exponentially with depth
@inline saturated_hydraulic_conductivity(hydraulics::SURFEXHydraulics, args...) = hydraulics.cond_sat

# TODO: Maybe we can borrow something better from SINDABD here; the SURFEX scheme is quite simplistic

@inline function mineral_porosity(hydraulics::SURFEXHydraulics, texture::SoilTexture)
    p₀ = hydraulics.porosity
    β_s = hydraulics.porosity_sand_coef
    por = p₀ + β_s*texture.sand*100
    return por
end

@inline function mineral_wilting_point(hydraulics::SURFEXHydraulics, texture::SoilTexture)
    β_w = hydraulics.wilting_point_coef
    wp = β_w*sqrt(texture.clay*100)
    return wp
end

@inline function mineral_field_capacity(hydraulics::SURFEXHydraulics, texture::SoilTexture)
    η = hydraulics.field_capacity_exp
    β_c = hydraulics.field_capacity_coef
    fc = β_c*(texture.clay*100)^η
    return fc
end
