"""
    BareGroundEvaporation{NF, GR} <: AbstractEvapotranspiration

Evaporation scheme for bare ground that calculates the humidity flux as

```math
E = \\beta \\frac{\\Delta q}{r_a}
```
where `Δq` is the vapor pressure deficit in terms of specific humidity, `rₐ` is aerodynamic resistance,
and `β` is an evaporation limiting factor.
"""
struct BareGroundEvaporation{NF, GR <: AbstractGroundEvaporationResistanceFactor} <: AbstractEvapotranspiration{NF}
    "Parameterization for ground resistance to evaporation/sublimation"
    ground_resistance::GR
end

BareGroundEvaporation(
    ::Type{NF};
    ground_resistance::GR = ConstantEvaporationResistanceFactor(one(NF))
) where {NF, GR} = BareGroundEvaporation{NF, GR}(; ground_resistance)

@propagate_inbounds surface_humidity_flux(i, j, grid, fields, evtr::BareGroundEvaporation) = fields.evaporation_ground[i, j]

# Process methods

variables(::BareGroundEvaporation) = (
    auxiliary(:evaporation_ground, XY(), units = u"m/s", desc = "Ground evaporation contribution to surface humidity flux"),
    input(:skin_temperature, XY(), units = u"°C", desc = "Skin temperature of the surface"),
)

function compute_auxiliary!(
        state, grid,
        evaporation::BareGroundEvaporation,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
        soil::Optional{AbstractSoil} = nothing
    )
    out = auxiliary_fields(state, evaporation)
    fields = get_fields(state, evaporation, atmos, soil; except = out)
    launch!(grid, XY, compute_evaporation_kernel!, out, fields, evaporation, atmos, constants, soil)
    return nothing
end

# Kernel functions

@propagate_inbounds function compute_evaporation_kernel!(
        out, i, j, grid, fields,
        evaporation::BareGroundEvaporation,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
        soil::Optional{AbstractSoil} = nothing
    )
    i, j = @index(Global, NTuple)
    Ts = fields.skin_temperature[i, j]
    rₐ = aerodynamic_resistance(i, j, grid, fields, aerodynamic_resistance) # aerodynamic resistance
    β = ground_evaporation_resistance_factor(i, j, grid, fields, evaporation.ground_resistance, soil)
    Δq = compute_humidity_vpd(i, j, grid, fields, atmos, constants, Ts)
    # Calculate water evaporation flux in m/s (positive upwards)
    return out.evaporation_ground[i, j, 1] = β * Δq / rₐ
end
