"""
    $TYPEDEF

Parameterization of soil porosity that simply specifies constant values
for the `mineral` and `organic` components.
"""
@parameterized @kwdef struct ConstantSoilPorosity{NF} <: AbstractSoilPorosity{NF}
    "Prescribed mineral soil porosity"
    @param mineral_porosity::NF = 0.49 (bounds = UnitInterval,)

    "Natural porosity of organic material"
    @param organic_porosity::NF = 0.9 (bounds = UnitInterval,)
end

ConstantSoilPorosity(::Type{NF}; kwargs...) where {NF} = ConstantSoilPorosity{NF}(; kwargs...)

@inline organic_porosity(params::ConstantSoilPorosity, texture::SoilTexture) = params.organic_porosity

@inline mineral_porosity(params::ConstantSoilPorosity, texture::SoilTexture) = params.mineral_porosity

"""
    $TYPEDEF

SURFEX parameterization of mineral soil porosity [noilhanISBA1996; Eq. (27)](@cite).

# References

* [noilhanISBA1996](@cite) Noilhan & Mahfouf, Global and Planetary Change (1996)
"""
@parameterized @kwdef struct SoilPorositySURFEX{NF} <: AbstractSoilPorosity{NF}
    "Assumed default porosity of the soil without sand"
    @param porosity_default::NF = 0.49 (bounds = UnitInterval,)

    "Linear coefficient of porosity adjustment for sand"
    @param porosity_sand_coef::NF = -0.11

    "Natural porosity of organic material"
    @param porosity_organic::NF = 0.9 (bounds = UnitInterval,)
end

SoilPorositySURFEX(::Type{NF}; kwargs...) where {NF} = SoilPorositySURFEX{NF}(; kwargs...)

@inline organic_porosity(params::SoilPorositySURFEX, texture::SoilTexture) = params.porosity_organic

@inline function mineral_porosity(params::SoilPorositySURFEX, texture::SoilTexture)
    p₀ = params.porosity_default
    β_s = params.porosity_sand_coef
    por = p₀ + β_s * texture.sand
    return por
end

# Kernel functions

function mineral_porosity(i, j, grid, fields, horizon::AbstractSoilHorizon)
    porosity_scheme = porosity(horizon)
    texture = soil_texture(i, j, grid, fields, horizon)
    return mineral_porosity(porosity_scheme, texture)
end

function organic_porosity(i, j, grid, fields, horizon::AbstractSoilHorizon)
    porosity_scheme = porosity(horizon)
    texture = soil_texture(i, j, grid, fields, horizon)
    return organic_porosity(porosity_scheme, texture)
end
