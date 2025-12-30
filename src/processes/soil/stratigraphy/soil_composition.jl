abstract type AbstractSoilMatrix{NF} end

organic_fraction(::AbstractSoilMatrix{NF}) where {NF} = zero(NF)

mineral_texture(solid::AbstractSoilMatrix) = solid.texture

"""
    $TYPEDEF

Represents the material composition of an elementary volume of soil.
The volume is decomposed into the key constitutents of water, ice, air,
and a mixture of organic and mineral solid material. 
"""
@kwdef struct SoilComposition{
    NF,
    T<:Union{NF, AbstractArray{NF}},
    Solid<:AbstractSoilMatrix
}
    "Natural porosity or void space of the soil"
    porosity::T = 0.5

    "Fraction of the soil pores occupied by water or ice"
    saturation::T = 1.0

    "Liquid (unfrozen) fraction of pore water"
    liquid::T = 1.0

    "Parmaeterization of the solid phase (matrix) of the soil"
    solid::Solid = MineralOrganic(texture = SoilTexture(), organic = 0.0)

    # Scalar constructor
    function SoilComposition(porosity::NF, saturation::NF, liquid::NF, solid::AbstractSoilMatrix{NF}) where {NF<:Number}
        @assert zero(NF) <= porosity <= one(NF)
        @assert zero(NF) <= saturation <= one(NF)
        @assert zero(NF) <= liquid <= one(NF)
        return new{NF, NF, typeof(solid)}(porosity, saturation, liquid, solid)
    end

    # Field/array constructor
    function SoilComposition(porosity::T, saturation::T, liquid::T, solid::AbstractSoilMatrix) where {NF, T<:AbstractArray{NF}}
        @assert size(porosity) == size(saturation) == size(liquid)
        return new{NF, T, typeof(solid)}(porosity, saturation, liquid, solid)
    end
end

@inline @propagate_inbounds function Base.getindex(soil::SoilComposition{<:AbstractArray}, idx...)
    return SoilComposition(
        soil.porosity[idx...],
        soil.saturation[idx...],
        soil.liquid[idx...],
        soil.solid[idx...]
    )
end

porosity(soil::SoilComposition) = soil.porosity

saturation(soil::SoilComposition) = soil.saturation

liquid_fraction(soil::SoilComposition) = soil.liquid

water_ice(soil::SoilComposition) = soil.porosity*soil.saturation

organic_fraction(soil::SoilComposition) = organic_fraction(soil.solid)

mineral_texture(soil::SoilComposition) = mineral_texture(soil.solid)

"""
    $SIGNATURES

Calculates the volumetric fractions of all constituents in the given soil volume
and returns them as a named tuple of the form `(; water, ice, air, solids...)`, where
`solids` correspodns to the volumetric fractions defined by the solid phase `soil.solid`.
"""
@inline function volumetric_fractions(soil::SoilComposition)
    # unpack relevant quantities
    let por = soil.porosity,
        sat = soil.saturation,
        liq = soil.liquid;
        # calculate volumetric fractions
        water_ice = sat * por
        water = water_ice * liq
        ice = water_ice * (1 - liq)
        air = (1 - sat) * por
        # get fractions of solid constituents
        solids = volumetric_fractions(soil.solid, solid_frac)
        return (; water, ice, air, solids...)
    end
end

"""
Soil matrix consisting of a simple, homogeneous mixture of mineral and organic material.
"""
struct MineralOrganic{T} <: AbstractSoilMatrix{T}
    "Mineral soil texture"
    texture::SoilTexture{T}

    "Organic soil fraction"
    organic::T
end

"""
Alias for `SoilComposition{T, MineralOrganic{T}}`
"""
const MineralOrganicSoil{T} = SoilComposition{T, MineralOrganic{T}}

@inline @propagate_inbounds function Base.getindex(solid::MineralOrganic{<:AbstractArray}, idx...)
    return MineralOrganic(
        solid.texture[idx...],
        solid.organic[idx...],
    )
end

@inline function volumetric_fractions(solid::MineralOrganic{FT}, solid_frac::FT) where {FT}
    organic = solid_frac * solid.organic
    mineral = solid_frac * (1 - solid.organic)
    return (; organic, mineral) 
end
