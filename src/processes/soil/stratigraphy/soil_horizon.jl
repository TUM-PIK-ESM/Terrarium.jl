"""
    $TYPEDEF

Represents an arbitrary soil horizon whose properties (texture and porosity) are assumed to
be constant across both space and time.
"""
@parameterized @kwdef struct ConstantSoilHorizon{NF, Porosity <: AbstractSoilPorosity{NF}} <: AbstractSoilHorizon{NF}
    "Material composition of mineral soil component"
    @component texture::SoilTexture{NF}

    "Parameterization of soil porosity"
    @component porosity::Porosity

    "Thickness of the soil horizon in meters"
    @param thickness::NF
end

function ConstantSoilHorizon(
        ::Type{NF};
        texture = SoilTexture(NF),
        porosity = ConstantSoilPorosity(NF),
        thickness = one(NF)
    ) where {NF}
    return ConstantSoilHorizon(texture, porosity, thickness)
end

@inline soil_texture(i, j, grid, fields, horizon::ConstantSoilHorizon) = horizon.texture

@inline soil_thickness(i, j, grid, fields, horizon::ConstantSoilHorizon) = horizon.thickness

"""
    $TYPEDEF

Represents an arbitrary soil horizon whose properties (texture and porosity) are prescribed
via input `Field`s and can therefore vary across space and (less commonly) time.
"""
@parameterized @kwdef struct PrescribedSoilHorizon{NF, Porosity <: AbstractSoilPorosity{NF}} <: AbstractSoilHorizon{NF}
    "Parameterization of soil porosity"
    @component porosity::Porosity
end

function PrescribedSoilHorizon(
        ::Type{NF};
        porosity = ConstantSoilPorosity(NF),
    ) where {NF}
    return PrescribedSoilHorizon(porosity)
end

variables(horizon::PrescribedSoilHorizon{NF}) where {NF} = (
    input(:sand, XY(), default = one(NF), bounds = UnitInterval, desc = "Mass fraction of sand in soil matrix"),
    input(:silt, XY(), bounds = UnitInterval, desc = "Mass fraction of silt in soil matrix"),
    input(:clay, XY(), bounds = UnitInterval, desc = "Mass fraction of clay in soil matrix"),
    input(:thickness, XY(), default = NF(Inf), bounds = Nonnegative, desc = "Thickness of soil horizon"),
)

@inline function soil_texture(i, j, grid, fields, horizon::PrescribedSoilHorizon{NF}) where {NF}
    sand = fields.sand[i, j]
    silt = fields.silt[i, j]
    clay = fields.clay[i, j]
    return SoilTexture(NF; sand, silt, clay)
end

@inline soil_thickness(i, j, grid, fields, horizon::PrescribedSoilHorizon{NF}) where {NF} = fields.thickness[i, j]
