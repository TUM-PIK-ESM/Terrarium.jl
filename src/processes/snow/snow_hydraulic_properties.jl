"""
    $TYPEDEF

Base type for snow hydraulic properties and parameterization schemes.
"""
abstract type AbstractSnowHydraulics{NF} end

"""
    $TYPEDEF

Constant, spatially homogeneous snow hydraulic properties: a saturated hydraulic conductivity `K_sat`
and a capillary retention `L_c`, setting the Darcy-type meltwater outflow (see
[`compute_meltwater_outflow`](@ref)). Default values follow [tarbotonSpatiallyDistributedEnergy1994](@cite).

Properties:
$TYPEDFIELDS

# References

* [tarbotonSpatiallyDistributedEnergy1994](@cite) Tarboton, Chowdhury and Jackson (1994)
"""
@parameterized @kwdef struct ConstantSnowHydraulics{NF} <: AbstractSnowHydraulics{NF}
    "Hydraulic conductivity at saturation"
    @param saturated_conductivity::NF = 20.0 / 3600.0 (units = u"m/s", bounds = Positive, scale = 1.0e-3)

    "Capillary retention `liq_c`: liquid fraction held against gravity before meltwater drains"
    @param capillary_retention::NF = 0.05 (bounds = UnitInterval,)
end

ConstantSnowHydraulics(::Type{NF}; kwargs...) where {NF} = ConstantSnowHydraulics{NF}(; kwargs...)
