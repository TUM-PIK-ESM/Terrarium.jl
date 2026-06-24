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

""" $TYPEDSIGNATURES """
@inline ground_evaporation_conductance(ET::AbstractEvapotranspiration, β, rₐ) = β / rₐ

"""
    $TYPEDSIGNATURES

Evaluate the bare-ground surface humidity flux [m/s] from the current skin temperature and
evaporation conductance in `fields`.
"""
@propagate_inbounds function compute_surface_humidity_flux(
        i, j, grid, fields,
        ::BareGroundEvaporation,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
    )
    Ts = fields.skin_temperature[i, j]
    g = fields.ground_evaporation_conductance[i, j]
    Δq = compute_specific_humidity_difference(i, j, grid, fields, atmos, constants, Ts)
    return g * Δq
end

# Top-level interface methods

variables(::BareGroundEvaporation) = (
    # Skin-driven vapor conductance β/rₐ (independent of skin temperature; held fixed during the SEB solve)
    auxiliary(:ground_evaporation_conductance, XY(), units = u"m/s", desc = "Ground evaporation vapor conductance"),
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

@propagate_inbounds function compute_evapotranspiration_conductances!(
        out, i, j, grid, fields,
        evaporation::BareGroundEvaporation,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil} = nothing
    )
    rₐ = aerodynamic_resistance(i, j, grid, fields, atmos) # aerodynamic resistance
    β = ground_evaporation_resistance_factor(i, j, grid, fields, evaporation.ground_resistance, soil)
    out.ground_evaporation_conductance[i, j, 1] = ground_evaporation_conductance(evaporation, β, rₐ)
    return out
end

@propagate_inbounds function compute_evapotranspiration_fluxes!(
        out, i, j, grid, fields,
        evaporation::BareGroundEvaporation,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
    )
    Ts = fields.skin_temperature[i, j]
    g = fields.ground_evaporation_conductance[i, j]
    # Evaporation flux at the current skin temperature. Re-running this kernel after the SEB solve
    # (see `LandModel`) refreshes it so it is consistent with the converged skin temperature.
    Δq = compute_specific_humidity_difference(i, j, grid, fields, atmos, constants, Ts)
    out.evaporation_ground[i, j, 1] = compute_evaporation_flux(evaporation, Δq, g)
    return out
end

# Kernels

@kernel inbounds = true function compute_auxiliary_kernel!(
        out, grid, fields,
        evapotranspiration::BareGroundEvaporation,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::AbstractSoil,
        args...
    )
    i, j = @index(Global, NTuple)

    # First compute conductances
    compute_evapotranspiration_conductances!(out, i, j, grid, fields, evapotranspiration, constants, atmos, soil, args...)
    # TODO: Annoyingly, we need to explicitly add these to `fields`; need a better solution to this problem
    conductances = (ground_evaporation_conductance = out.ground_evaporation_conductance,)
    fields = merge(fields, conductances)
    # Compute ET fluxes from stored conductances
    compute_evapotranspiration_fluxes!(out, i, j, grid, fields, evapotranspiration, constants, atmos)
end
