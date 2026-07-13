"""
    $TYPEDEF

Base type for snow hydraulic properties and parameterization schemes.
"""
abstract type AbstractSnowHydraulics{NF} end

@parameterized @kwdef struct ConstantSnowHydraulics{NF} <: AbstractSnowHydraulics{NF}
    "Hydraulic conductivity at saturation"
    @param saturated_conductivity::NF = 1.0e-5 (units = u"m/s", bounds = Positive, scale = 1.0e-5)

    "Capillary retention `L_c`: liquid fraction held against gravity before meltwater drains"
    @param capillary_retention::NF = 0.05 (bounds = UnitInterval,)
end

ConstantSnowHydraulics(::Type{NF}; kwargs...) where {NF} = ConstantSnowHydraulics{NF}(; kwargs...)
