"""
    $TYPEDEF

Represents an arbitrary soil horizon whose properties (texture and porosity) are assumed to
be constant across both space and time.
"""
@parameterized struct ConstantSoilHorizon{NF, name, Porosity <: AbstractSoilPorosity{NF}} <: AbstractSoilHorizon{NF, name}
    "Material composition of mineral soil component"
    @component texture::SoilTexture{NF}

    "Parameterization of soil porosity"
    @component porosity::Porosity

    "Thickness of the soil horizon in meters"
    @param thickness::NF

    function ConstantSoilHorizon(name::Symbol, texture::SoilTexture, porosity::AbstractSoilPorosity, thickness::NF) where {NF}
        return new{NF, name, typeof(porosity)}(texture, porosity, thickness)
    end
end

function ConstantSoilHorizon(
        ::Type{NF}, name::Symbol;
        texture = SoilTexture(NF),
        porosity = ConstantSoilPorosity(NF),
        thickness = one(NF)
    ) where {NF}
    return ConstantSoilHorizon(name, texture, porosity, thickness)
end

@inline soil_texture(i, j, grid, fields, horizon::ConstantSoilHorizon) = horizon.texture

@inline soil_thickness(i, j, grid, fields, horizon::ConstantSoilHorizon) = horizon.thickness

"""
    $TYPEDEF

Represents an arbitrary soil horizon whose properties (texture and porosity) are prescribed
via input `Field`s and can therefore vary across space and (less commonly) time.
"""
@parameterized struct PrescribedSoilHorizon{NF, name, Porosity <: AbstractSoilPorosity{NF}} <: AbstractSoilHorizon{NF, name}
    "Parameterization of soil porosity"
    @component porosity::Porosity

    function PrescribedSoilHorizon(name::Symbol, porosity::AbstractSoilPorosity{NF}) where {NF}
        return new{NF, name, typeof(porosity)}(porosity)
    end
end

function PrescribedSoilHorizon(
        ::Type{NF}, name::Symbol;
        porosity = ConstantSoilPorosity(NF),
    ) where {NF}
    return PrescribedSoilHorizon(name, porosity)
end

variables(::PrescribedSoilHorizon{NF}) where {NF} = (
    input(:sand_fraction, XY(), default = one(NF), bounds = UnitInterval, desc = "Mass fraction of sand in soil matrix"),
    input(:silt_fraction, XY(), bounds = UnitInterval, desc = "Mass fraction of silt in soil matrix"),
    input(:clay_fraction, XY(), bounds = UnitInterval, desc = "Mass fraction of clay in soil matrix"),
    input(:thickness, XY(), default = NF(Inf), bounds = Nonnegative, desc = "Thickness of soil horizon"),
)

@inline function soil_texture(i, j, grid, fields, horizon::PrescribedSoilHorizon{NF}) where {NF}
    sand = fields.sand_fraction[i, j]
    silt = fields.silt_fraction[i, j]
    clay = fields.clay_fraction[i, j]
    return SoilTexture(NF; sand, silt, clay)
end

@inline soil_thickness(i, j, grid, fields, horizon::PrescribedSoilHorizon{NF}) where {NF} = fields.thickness[i, j]
