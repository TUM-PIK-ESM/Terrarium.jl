"""
    BareGroundEvaporation{NF, GR} <: AbstractEvapotranspiration

Evaporation scheme for bare ground that calculates the humidity flux as

```math
E = \\beta \\frac{\\Delta q}{r_a}
```
where `Δq` is the specific humidity difference, `rₐ` is aerodynamic resistance,
and `β` is an evaporation limiting factor.
"""
struct BareGroundEvaporation{NF, GR <: AbstractGroundEvaporationResistanceFactor} <: AbstractEvapotranspiration{NF}
    "Parameterization for ground resistance to evaporation/sublimation"
    ground_resistance::GR
end

BareGroundEvaporation(
    ::Type{NF};
    ground_resistance::GR = SoilMoistureResistanceFactor(NF)
) where {NF, GR} = BareGroundEvaporation{NF, GR}(ground_resistance)

@propagate_inbounds surface_humidity_flux(i, j, grid, fields, evaporation::BareGroundEvaporation, args...) = fields.evaporation_ground[i, j]

"""
    $TYPEDSIGNATURES

Evaluate the bare-ground surface humidity flux [m/s] from the skin temperature read live from
`fields`, reusing the precomputed evaporation conductance `β/rₐ`. The whole flux is skin-driven,
so this is simply `(β/rₐ)·Δq(Tₛ)`. Used by the SEB during the implicit skin-temperature solve so
that the latent heat flux responds to the skin temperature.
"""
@propagate_inbounds function surface_humidity_flux(
        i, j, grid, fields,
        ::BareGroundEvaporation,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
    )
    Ts = fields.skin_temperature[i, j]
    g = fields.evaporation_conductance[i, j]
    Δq = compute_specific_humidity_difference(i, j, grid, fields, atmos, constants, Ts)
    return g * Δq
end

# Top-level interface methods

variables(::BareGroundEvaporation) = (
    # Skin-driven vapor conductance β/rₐ (independent of skin temperature; held fixed during the SEB solve)
    auxiliary(:evaporation_conductance, XY(), units = u"m/s", desc = "Ground evaporation vapor conductance"),
    auxiliary(:evaporation_ground, XY(), units = u"m/s", desc = "Ground evaporation contribution to surface humidity flux"),
    input(:skin_temperature, XY(), units = u"°C", desc = "Skin temperature of the surface"),
)

""" $TYPEDSIGNATURES """
function compute_auxiliary!(
        state, grid,
        evaporation::BareGroundEvaporation,
        ::NoCanopyInterception,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil} = nothing,
        args...
    )
    out = auxiliary_fields(state, evaporation)
    fields = get_fields(state, evaporation, atmos, soil; except = out)
    launch!(grid, XY, compute_auxiliary_kernel!, out, fields, evaporation, constants, atmos, soil)
    return nothing
end

# Kernel functions

@propagate_inbounds function compute_evapotranspiration!(
        out, i, j, grid, fields,
        evaporation::BareGroundEvaporation,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil} = nothing
    )
    Ts = fields.skin_temperature[i, j]
    rₐ = aerodynamic_resistance(i, j, grid, fields, atmos) # aerodynamic resistance
    β = ground_evaporation_resistance_factor(i, j, grid, fields, evaporation.ground_resistance, soil)
    # Store the skin-driven evaporation conductance β/rₐ (independent of skin temperature). This is
    # the source of truth from which the SEB derives the latent heat flux lazily during the skin
    # temperature solve (see the lazy `surface_humidity_flux` method above).
    g = β / rₐ
    out.evaporation_conductance[i, j, 1] = g
    # Evaporation flux at the current skin temperature. Re-running this kernel after the SEB solve
    # (see `LandModel`) refreshes it so it is consistent with the converged skin temperature.
    Δq = compute_specific_humidity_difference(i, j, grid, fields, atmos, constants, Ts)
    out.evaporation_ground[i, j, 1] = g * Δq
    return out
end
