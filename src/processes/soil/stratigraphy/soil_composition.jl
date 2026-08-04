"""
    $TYPEDEF

Represents the material composition of an elementary volume of soil.
The volume is decomposed into the key constitutents of water, ice, air,
and a mixture of organic and mineral solid material.

Properties:
$FIELDS
"""
struct SoilComposition{NF, Solid <: AbstractSoilMatrix{NF}}
    "Natural porosity or void space of the soil"
    porosity::NF

    "Fraction of the soil pores occupied by water or ice"
    saturation::NF

    "Liquid (unfrozen) fraction of pore water"
    liquid::NF

    "Parameterization of the solid phase (matrix) of the soil"
    solid::Solid

    function SoilComposition(porosity::NF, saturation::NF, liquid::NF, solid::AbstractSoilMatrix{NF}) where {NF <: Number}
        return new{NF, typeof(solid)}(porosity, saturation, liquid, solid)
    end
end

"""
    $SIGNATURES

Validating keyword constructor for `SoilVolume`; enforces that `porosity`, `saturation`,
and `liquid` lie in the unit interval.
"""
function SoilComposition(;
        porosity = 0.5,
        saturation = 1.0,
        liquid = 1.0,
        solid = MineralOrganic(texture = SoilTexture(), organic = 0.0)
    )
    @assert zero(porosity) <= porosity <= one(porosity)
    @assert zero(saturation) <= saturation <= one(saturation)
    @assert zero(liquid) <= liquid <= one(liquid)
    return SoilComposition(porosity, saturation, liquid, solid)
end

@inline porosity(soil::SoilComposition) = soil.porosity

@inline saturation(soil::SoilComposition) = soil.saturation

@inline liquid_fraction(soil::SoilComposition) = soil.liquid

@inline water_ice(soil::SoilComposition) = soil.porosity * soil.saturation

@inline water(soil::SoilComposition) = soil.liquid * soil.porosity * soil.saturation

@inline organic_fraction(soil::SoilComposition) = organic_fraction(soil.solid)

@inline mineral_texture(soil::SoilComposition) = mineral_texture(soil.solid)

"""
    $TYPEDSIGNATURES

Calculates the volumetric fractions of all constituents in the given soil volume
and returns them as a named tuple of the form `(; water, ice, air, solids...)`, where
`solids` corresponds to the volumetric fractions defined by the solid phase `soil.solid`.
"""
@inline function volumetric_fractions(soil::SoilComposition)
    # unpack relevant quantities
    let por = soil.porosity,
            sat = soil.saturation,
            liq = soil.liquid
        # calculate volumetric fractions
        water_ice = sat * por
        water = water_ice * liq
        ice = water_ice * (1 - liq)
        air = (1 - sat) * por
        solid_frac = 1 - por
        # get fractions of solid constituents
        solids = volumetric_fractions(soil.solid, solid_frac)
        return (; water, ice, air, solids...)
    end
end

"""
    $TYPEDEF

Soil matrix consisting of a simple, homogeneous mixture of mineral and organic material.

Properties:
$TYPEDFIELDS
"""
struct MineralOrganic{NF} <: AbstractSoilMatrix{NF}
    "Mineral soil texture"
    texture::SoilTexture{NF}

    "Organic soil fraction"
    organic::NF

    function MineralOrganic(texture::SoilTexture{NF}, organic::NF) where {NF}
        return new{NF}(texture, organic)
    end
end

"""
    $SIGNATURES

Validating keyword constructor for `MineralOrganic`; enforces that `organic` lies in the
unit interval.
"""
function MineralOrganic(; texture = SoilTexture(), organic = zero(eltype(texture)))
    @assert zero(organic) <= organic <= one(organic) "organic content must be between zero and one"
    return MineralOrganic(texture, organic)
end

"""
Alias for `SoilComposition{T, MineralOrganic{T}}`
"""
const MineralOrganicSoil{NF} = SoilComposition{NF, MineralOrganic{NF}}

@inline mineral_texture(solid::MineralOrganic) = solid.texture

@inline organic_fraction(solid::MineralOrganic) = solid.organic

"""
    $TYPEDSIGNATURES

Compute the volumetric fractions of the solid phase scaled by the overall solid fraction
of the soil `solid_frac`.
"""
@inline function volumetric_fractions(solid::MineralOrganic{NF}, solid_frac::NF) where {NF}
    organic = solid_frac * solid.organic
    mineral = solid_frac * (1 - solid.organic)
    return (; organic, mineral)
end
