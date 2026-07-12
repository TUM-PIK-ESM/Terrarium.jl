"""
    $TYPEDEF

Constant, spatially homogeneous bulk snow density `ρ_s`. This is the default (and currently only) snow
density scheme for [`SingleLayerSnow`](@ref).

Properties:
$TYPEDFIELDS
"""
@parameterized @kwdef struct ConstantSnowDensity{NF} <: AbstractSnowDensity
    "Bulk snow density `ρ_s`"
    @param density::NF = 250.0 (units = u"kg/m^3", bounds = Positive)
end

ConstantSnowDensity(::Type{NF}; kwargs...) where {NF} = ConstantSnowDensity{NF}(; kwargs...)

"""
    $TYPEDSIGNATURES

Return the constant bulk snow density `ρ_s` [kg/m³].
"""
@inline snow_density(density::ConstantSnowDensity) = density.density

@propagate_inbounds compute_snow_density(i, j, grid, fields, density::ConstantSnowDensity) = snow_density(density)
