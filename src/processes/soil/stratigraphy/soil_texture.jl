"""
    $TYPEDEF

Represents soil texture as a fractional mixture of sand, silt, and clay. Accepts values 
"""
@kwdef struct SoilTexture{NF} <: AbstractSoilTexture{NF}
    sand::NF = 1.0
    clay::NF = 0.0
    silt::NF = 1 - sand - clay

    SoilTexture(sand, silt, clay) = SoilTexture{promote_type(typeof(sand), typeof(silt), typeof(clay))}(sand, silt, clay)
    function SoilTexture{NF}(sand, silt, clay) where {NF<:Number}
        @assert zero(sand) <= sand <= one(sand)
        @assert zero(silt) <= silt <= one(silt)
        @assert zero(clay) <= clay <= one(clay)
        @assert sand + silt + clay ≈ 1.0 "sand, silt, and clay fractions must sum to unity"
        return new{NF}(sand, silt, clay)
    end
end

SoilTexture(::Type{NF}; kwargs...) where {NF} = SoilTexture{NF}(; kwargs...)

Base.eltype(::SoilTexture{NF}) where {NF} = NF

@inline @propagate_inbounds function Base.getindex(texture::SoilTexture{<:AbstractArray}, idx...)
    return SoilTexture(
        texture.sand[idx...],
        texture.silt[idx...],
        texture.clay[idx...]
    )
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
SoilTexture(name::Symbol) = SoilTexture(Float64, Val{name}())
SoilTexture(::Type{NF}, name::Symbol) where {NF} = SoilTexture(NF, Val{name}())
SoilTexture(::Type{NF}, ::Val{:sand}) where {NF} = SoilTexture(NF, sand=1.0, silt=0.0, clay=0.0)
SoilTexture(::Type{NF}, ::Val{:silt}) where {NF} = SoilTexture(NF, sand=0.0, silt=1.0, clay=0.0)
SoilTexture(::Type{NF}, ::Val{:clay}) where {NF} = SoilTexture(NF, sand=0.0, silt=0.0, clay=1.0)
SoilTexture(::Type{NF}, ::Val{:sandyclay}) where {NF} = SoilTexture(NF, sand=0.50, silt=0.0, clay=0.50)
SoilTexture(::Type{NF}, ::Val{:siltyclay}) where {NF} = SoilTexture(NF, sand=0.0, silt=0.50, clay=0.50)
SoilTexture(::Type{NF}, ::Val{:loam}) where {NF} = SoilTexture(NF, sand=0.40, silt=0.40, clay=0.20)
SoilTexture(::Type{NF}, ::Val{:sandyloam}) where {NF} = SoilTexture(NF, sand=0.80, silt=0.10, clay=0.10)
SoilTexture(::Type{NF}, ::Val{:siltyloam}) where {NF} = SoilTexture(NF, sand=0.10, silt=0.80, clay=0.10)
SoilTexture(::Type{NF}, ::Val{:clayloam}) where {NF} = SoilTexture(NF, sand=0.30, silt=0.30, clay=0.40)
