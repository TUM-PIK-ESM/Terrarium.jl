"""
    SoilTexture{NF}

Represents soil texture as a fractional mixture of sand, silt, and clay.
"""
@kwdef struct SoilTexture{NF}
    sand::NF = 1.0
    clay::NF = 0.0
    silt::NF = 1 - sand - clay
    function SoilTexture(sand::NF, silt::NF, clay::NF) where {NF}
        @assert zero(sand) <= sand <= one(sand)
        @assert zero(silt) <= silt <= one(silt)
        @assert zero(clay) <= clay <= one(clay)
        @assert sand + silt + clay ≈ 1.0 "sand, silt, and clay fractions must sum to unity"
        return new{NF}(sand, silt, clay)
    end
end

function Base.convert(::SoilTexture{NewType}, texture::SoilTexture) where {NewType}
    return SoilTexture(
        convert(NewType, texture.sand),
        convert(NewType, texture.silt),
        convert(NewType, texture.clay)
    )
end

# Presets for common soil textures
# Borrowed from CryoGrid.jl:
# https://github.com/CryoGrid/CryoGrid.jl/blob/master/src/Physics/Soils/soil_texture.jl
SoilTexture(name::Symbol) = SoilTexture(Val{name}())
SoilTexture(::Val{:sand}) = SoilTexture(sand=1.0, silt=0.0, clay=0.0)
SoilTexture(::Val{:silt}) = SoilTexture(sand=0.0, silt=1.0, clay=0.0)
SoilTexture(::Val{:clay}) = SoilTexture(sand=0.0, silt=0.0, clay=1.0)
SoilTexture(::Val{:sandyclay}) = SoilTexture(sand=0.50, silt=0.0, clay=0.50)
SoilTexture(::Val{:siltyclay}) = SoilTexture(sand=0.0, silt=0.50, clay=0.50)
SoilTexture(::Val{:loam}) = SoilTexture(sand=0.40, silt=0.40, clay=0.20)
SoilTexture(::Val{:sandyloam}) = SoilTexture(sand=0.80, silt=0.10, clay=0.10)
SoilTexture(::Val{:siltyloam}) = SoilTexture(sand=0.10, silt=0.80, clay=0.10)
SoilTexture(::Val{:clayloam}) = SoilTexture(sand=0.30, silt=0.30, clay=0.40)
