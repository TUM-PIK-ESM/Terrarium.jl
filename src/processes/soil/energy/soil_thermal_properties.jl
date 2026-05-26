# Note: Should sand, silt, and clay have separate thermal properties?
"""
    $TYPEDEF

Properties:
$TYPEDFIELDS

Default values from [hillelIntroductionSoilPhysics1982](@cite).

# References

* [hillelIntroductionSoilPhysics1982](@cite) Hillel, Academic Press (1982)
"""
@parameterized @kwdef struct SoilThermalConductivities{NF}
    "Thermal conductivity of water"
    @param water::NF = 0.57 (units = u"W/m/K", bounds = Positive)

    "Thermal conductivity of ice"
    @param ice::NF = 2.2 (units = u"W/m/K", bounds = Positive)

    "Thermal conductivity of air"
    @param air::NF = 0.025 (units = u"W/m/K", bounds = Positive)

    "Thermal conductivity of mineral soil constituents"
    @param mineral::NF = 3.8 (units = u"W/m/K", bounds = Positive)

    "Thermal conductivity of organic soil constituents"
    @param organic::NF = 0.25 (units = u"W/m/K", bounds = Positive)
end

SoilThermalConductivities(::Type{NF}; kwargs...) where {NF} = SoilThermalConductivities{NF}(; kwargs...)

"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
@parameterized @kwdef struct SoilHeatCapacities{NF}
    "Volumetric heat capacity of water"
    @param water::NF = 4.2e6 (units = u"J/K/m^3", bounds = Positive, scale = 1.0e6)

    "Volumetric heat capacity of ice"
    @param ice::NF = 1.9e6 (units = u"J/K/m^3", bounds = Positive, scale = 1.0e6)

    "Volumetric heat capacity of air"
    @param air::NF = 0.00125e6 (units = u"J/K/m^3", bounds = Positive, scale = 1.0e6)

    "Volumetric heat capacity of mineral soil"
    @param mineral::NF = 2.0e6 (units = u"J/K/m^3", bounds = Positive, scale = 1.0e6)

    "Volumetric heat capacity of organic soil"
    @param organic::NF = 2.5e6 (units = u"J/K/m^3", bounds = Positive, scale = 1.0e6)
end

SoilHeatCapacities(::Type{NF}; kwargs...) where {NF} = SoilHeatCapacities{NF}(; kwargs...)

# TODO: In principle, these types could change for different soil parameterizations.
# This is something we should ideally allow for.
"""
    $TYPEDEF

Properties:
$TYPEDFIELDS
"""
@parameterized @kwdef struct SoilThermalProperties{NF, FC, CondBulk}
    "Thermal conductivities for all constituents"
    @component conductivities::SoilThermalConductivities{NF}

    "Method for computing bulk thermal conductivity from constituents"
    @component bulk_conductivity::CondBulk

    "Thermal conductivities for all constituents"
    @component heat_capacities::SoilHeatCapacities{NF}

    "Freezing characteristic curve needed for energy-temperature closure"
    @component freezecurve::FC
end

SoilThermalProperties(
    ::Type{NF};
    conductivities::SoilThermalConductivities{NF} = SoilThermalConductivities(NF),
    bulk_conductivity::AbstractBulkWeighting = InverseQuadratic(),
    heat_capacities::SoilHeatCapacities{NF} = SoilHeatCapacities(NF),
    freezecurve::FreezeCurve = FreeWater()
) where {NF} = SoilThermalProperties{NF, typeof(freezecurve), typeof(bulk_conductivity)}(conductivities, bulk_conductivity, heat_capacities, freezecurve)

freezecurve(
    ::SoilThermalProperties{NF, FreeWater},
    ::AbstractSoilHydrology
) where {NF} = FreeWater()

"""
    $SIGNATURES

Compute the bulk thermal conductivity of the given soil volume.
"""
@inline function compute_thermal_conductivity(props::SoilThermalProperties, soil::SoilVolume)
    κs = getproperties(props.conductivities)
    fracs = volumetric_fractions(soil)
    # apply bulk conductivity weighting
    return props.bulk_conductivity(κs, fracs)
end

"""
    $SIGNATURES

Compute the bulk heat capacity of the given soil volume.
"""
@inline function heat_capacity(props::SoilThermalProperties, soil::SoilVolume)
    cs = getproperties(props.heat_capacities)
    fracs = volumetric_fractions(soil)
    # for heat capacity, we just do a weighted average
    return sum(fastmap(*, cs, fracs))
end

"""
The inverse quadratic (or "quadratic parallel") bulk thermal conductivity formula ([cosenzaSimultaneousDeterminationThermal2003](@cite)):

```math
k = \\left[\\sum_{i=1}^N θᵢ\\sqrt{kᵢ}\\right]^2
```

# References
* [cosenzaSimultaneousDeterminationThermal2003](@cite) Cosenza et al., European Journal of Soil Science (2003)
"""
struct InverseQuadratic <: AbstractBulkWeighting end

(f::InverseQuadratic)(x::Real, weight::Real) = sqrt(x) * weight
# we use fastmap here so that the ordering of named tuples can be arbitrary
(f::InverseQuadratic)(xs, weights) = sum(fastmap(f, xs, weights))^2
