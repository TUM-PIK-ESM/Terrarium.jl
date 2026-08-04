"""
    $TYPEDEF

Constant, spatially homogeneous bulk snow density `ρ_snow`. This is the default (and currently only) snow
density scheme for [`SingleLayerSnow`](@ref). Default bulk snow density follows
[westermannCryoGrid3Simulating2016](@cite).

Properties:
$TYPEDFIELDS

# References

* [westermannCryoGrid3Simulating2016](@cite) Westermann et al., Geoscientific Model Development (2016)
"""
@parameterized @kwdef struct ConstantSnowDensity{NF} <: AbstractSnowDensity{NF}
    "Bulk snow density `ρ_snow`"
    @param density::NF = 250.0 (units = u"kg/m^3", bounds = Positive)
end

ConstantSnowDensity(::Type{NF}; kwargs...) where {NF} = ConstantSnowDensity{NF}(; kwargs...)

"""
    $TYPEDSIGNATURES

Return the constant bulk snow density `ρ_snow` [kg/m³].
"""
@inline snow_density(density::ConstantSnowDensity) = density.density

@propagate_inbounds compute_snow_density(i, j, grid, fields, density::ConstantSnowDensity) = snow_density(density)
