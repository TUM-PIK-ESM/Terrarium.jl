abstract type AbstractSoilComposition end

"""
    $TYPEDEF

Represents the material composition of an elementary volume of soil.
The volume is decomposed into the key constitutents of water, ice, air,
and a mixture of organic and mineral solid material. 
"""
@kwdef struct SoilComposition{NF} <: AbstractSoilComposition
    "Natural porosity or void space of the soil"
    porosity::NF = 0.5

    "Fraction of the soil pores occupied by water or ice"
    saturation::NF = 1.0

    "Liquid (unfrozen) fraction of pore water"
    liquid::NF = 1.0

    "Fraction of the soil "
    organic::NF = 0.0

    "Soil texture (mineral solid constituents)"
    texture::SoilTexture{NF} = SoilTexture()
end

porosity(soil::SoilComposition) = soil.porosity

saturation(soil::SoilComposition) = soil.saturation

liquid_fraction(soil::SoilComposition) = soil.liquid

organic_fraction(soil::SoilComposition) = soil.organic

"""
    volumetric_fractions(soil::SoilComposition)

Calculates the volumetric fractions of all constituents in the given soil volume
and returns them as a named tuple of the form `(; water, ice, air, mineral, organic)`.
"""
@inline function volumetric_fractions(soil::SoilComposition)
    # unpack relevant quantities
    let por = soil.porosity,
        sat = soil.saturation,
        liq = soil.liquid,
        org = soil.organic;
        # calculate volumetric fractions
        water_ice = sat*por
        water = water_ice*liq
        ice = water_ice*(1-liq)
        air = (1-sat)*por
        mineral = (1-por)*(1-org)
        organic = (1-por)*org
        return (; water, ice, air, mineral, organic)
    end
end
