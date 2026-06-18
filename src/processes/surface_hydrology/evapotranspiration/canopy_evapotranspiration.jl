"""
    $TYPEDEF

Canopy evapotranspiration scheme from PALADYN ([willeitPALADYNV10Comprehensive2016; Eq. (5)](@cite)) 
that includes a canopy evaporation term based on the saturation fraction of canopy water defined by the
canopy hydrology scheme.

```math
E_{\\text{ground}} = \\beta \\frac{\\Delta q}{r_a + r_e}
```
```math
E_{\\text{can}} = f_{\\text{can}} \\frac{\\Delta q}{r_a}
```
```math
T_{\\text{can}} = \\frac{\\Delta q}{r_a + r_s}
```

Properties:
$FIELDS

# References
* [willeitPALADYNV10Comprehensive2016](@cite) Willeit and Ganopolski, Geoscientific Model Development (2016)
"""
@parameterized @kwdef struct PALADYNCanopyEvapotranspiration{NF, GR <: AbstractGroundEvaporationResistanceFactor} <: AbstractEvapotranspiration{NF}
    "Drag coefficient for the transfer of heat and water between the ground and canopy"
    @param C_can::NF (bounds = Positive,)

    "Parameterization for ground resistance to evaporation/sublimation"
    @component ground_resistance::GR
end

function PALADYNCanopyEvapotranspiration(
        ::Type{NF};
        C_can = NF(0.006),
        ground_resistance = ConstantEvaporationResistanceFactor(typeof(C_can))
    ) where {NF}
    return PALADYNCanopyEvapotranspiration{NF, typeof(ground_resistance)}(C_can, ground_resistance)
end

# TODO: The following ET functions all have the same basic functional form and can be generalized to a function
# that takes the humidity gradient/difference and an arbitrary number of resistance terms. We could consider
# making such an abstraction, but we should first consider what exactly the benefits would be since it also
# obfuscates that actual calculations, and the equations are quite simple. Perhaps one benefit would be reduced
# unit testing overhead?

"""
    $TYPEDSIGNATURES

Compute the transpiration vapor conductance [m/s] from aerodynamic resistance `rₐ` and stomatal
conductance `gw_can`. The transpiration flux is this conductance times the humidity gradient.
"""
@inline function transpiration_conductance(::PALADYNCanopyEvapotranspiration{NF}, rₐ, gw_can) where {NF}
    rₛ = 1 / max(gw_can, sqrt(eps(NF)))  # stomatal resistance as reciprocal of conductance
    return 1 / (rₐ + rₛ)
end

"""
    $TYPEDSIGNATURES

Compute the canopy evaporation vapor conductance [m/s] from canopy saturation fraction `f_can`
and aerodynamic resistance `rₐ`. The canopy evaporation flux is this conductance times the
humidity gradient.
"""
@inline canopy_evaporation_conductance(::PALADYNCanopyEvapotranspiration, f_can, rₐ) = f_can / rₐ

"""
    $TYPEDSIGNATURES

Compute transpiration from the given humidity gradient, aerodynamic resistance `rₐ` and stomatal conductance `gw_can`.
"""
@inline function compute_transpiration(evtr::PALADYNCanopyEvapotranspiration, Δq, rₐ, gw_can)
    # Calculate transpiration flux in m/s (positive upwards)
    E_trp = transpiration_conductance(evtr, rₐ, gw_can) * Δq
    return E_trp
end

"""
    $TYPEDSIGNATURES

Compute evaporation from the ground below the canopy, following [willeitPALADYNV10Comprehensive2016; Eq. (5)](@cite);
`Δq` is the humidity gradient, `β` is the ground evaporation resistance factor, `rₐ` is aerodynamic resistance,
and `rₑ` is aerodynamic resistance between the ground and canopy.

# References
* [willeitPALADYNV10Comprehensive2016](@cite) Willeit and Ganopolski, Geoscientific Model Development (2016)
"""
@inline function compute_evaporation_ground(::PALADYNCanopyEvapotranspiration, Δq, β, rₐ, rₑ)
    # Calculate ground evaporation flux in m/s (positive upwards)
    E_gnd = β * Δq / (rₐ + rₑ)
    return E_gnd
end

"""
    $TYPEDSIGNATURES

Compute evaporation of water intercepted by the canopy from humidity gradient `Δq`, canopy saturation fraction
`f_can`, and aerodynamic resistance `rₐ`.
"""
@inline function compute_evaporation_canopy(evtr::PALADYNCanopyEvapotranspiration, Δq, f_can, rₐ)
    # Calculate canopy evaporation flux in m/s (positive upwards)
    E_can = canopy_evaporation_conductance(evtr, f_can, rₐ) * Δq
    return E_can
end

# Top-level interface methods

variables(::PALADYNCanopyEvapotranspiration{NF}) where {NF} = (
    # Skin-driven vapor conductances (independent of skin temperature; held fixed during the SEB solve)
    auxiliary(:canopy_evaporation_conductance, XY(); units = u"m/s", desc = "Canopy evaporation vapor conductance"),
    auxiliary(:transpiration_conductance, XY(); units = u"m/s", desc = "Transpiration vapor conductance"),
    # Partitioned humidity fluxes (skin-driven terms refreshed from the converged skin temperature by the finalize pass)
    auxiliary(:evaporation_canopy, XY(); desc = "Canopy evaporation contribution to surface humidity flux", units = u"m/s"),
    auxiliary(:evaporation_ground, XY(), units = u"m/s", desc = "Ground evaporation contribution to surface humidity flux"),
    auxiliary(:transpiration, XY(), units = u"m/s", desc = "Transpiration contribution to surface humidity flux"),
    input(:skin_temperature, XY(); units = u"°C", desc = "Skin temperature"),
    input(:ground_temperature, XY(); default = NF(1), units = u"°C", desc = "Ground surface temperature"),
)

@propagate_inbounds function surface_humidity_flux(i, j, grid, fields, ::PALADYNCanopyEvapotranspiration)
    E_gnd = fields.evaporation_ground[i, j]
    E_can = fields.evaporation_canopy[i, j]
    T_can = fields.transpiration[i, j]
    return E_gnd + E_can + T_can
end

"""
    $TYPEDSIGNATURES

Evaluate the surface humidity flux [m/s] from the skin temperature read live from `fields`,
reusing the precomputed vapor conductances. The skin-driven terms (canopy evaporation and
transpiration) are recomputed from `Δq(Tₛ)`; the ground evaporation term is held fixed since it
is driven by the ground temperature, not the skin temperature. This method is used by the SEB
during the implicit skin-temperature solve, where `fields.skin_temperature` is the current
iterate, so the latent heat flux responds to the skin temperature.
"""
@propagate_inbounds function surface_humidity_flux(
        i, j, grid, fields,
        ::PALADYNCanopyEvapotranspiration,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
    )
    Ts = fields.skin_temperature[i, j]
    g_can = fields.canopy_evaporation_conductance[i, j]
    g_trp = fields.transpiration_conductance[i, j]
    E_gnd = fields.evaporation_ground[i, j]
    Δqs = compute_specific_humidity_difference(i, j, grid, fields, atmos, constants, Ts)
    return (g_can + g_trp) * Δqs + E_gnd
end

""" $TYPEDSIGNATURES """
function compute_auxiliary!(
        state, grid,
        evapotranspiration::PALADYNCanopyEvapotranspiration,
        canopy_interception::AbstractCanopyInterception,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::AbstractSoil,
        vegetation::AbstractVegetation,
        args...
    )
    out = auxiliary_fields(state, evapotranspiration)
    fields = get_fields(state, evapotranspiration, canopy_interception, atmos, soil, vegetation; except = out)
    launch!(grid, XY, compute_auxiliary_kernel!, out, fields, evapotranspiration, canopy_interception, constants, atmos, soil, vegetation)
    return nothing
end

# Kernel functions

"""
    $TYPEDEF

Compute `transpiration`, `evaporation_ground`, and `evaporation_canopy` fluxes on `grid`
for the given scheme `evapotranspiration` and process dependencies.
"""
@propagate_inbounds function compute_evapotranspiration!(
        out, i, j, grid, fields,
        evapotranspiration::PALADYNCanopyEvapotranspiration,
        canopy_interception::AbstractCanopyInterception,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::AbstractSoil,
        vegetation::AbstractVegetation,
        args...
    )
    # Get inputs
    Ts = fields.skin_temperature[i, j] # skin temperature (top of canopy)
    Tg = fields.ground_temperature[i, j] # ground temperature (top snow/soil layer)
    gw_can = fields.canopy_water_conductance[i, j] # stomatal conductance (assumed to be defined by vegetation)

    # Compute VPD and resistance terms
    Δqs = compute_specific_humidity_difference(i, j, grid, fields, atmos, constants, Ts) # humidity difference between canopy and atmosphere
    Δqg = compute_specific_humidity_difference(i, j, grid, fields, atmos, constants, Tg) # humidity difference between ground and canopy
    rₐ = aerodynamic_resistance(i, j, grid, fields, atmos) # aerodynamic resistance
    rₑ = aerodynamic_resistance(i, j, grid, fields, atmos, evapotranspiration, vegetation) # aerodynamic resistance between ground and canopy
    f_can = saturation_canopy_water(i, j, grid, fields, canopy_interception)
    β = ground_evaporation_resistance_factor(i, j, grid, fields, evapotranspiration.ground_resistance, soil)

    # Store skin-driven vapor conductances (independent of skin temperature). These are the
    # source of truth from which the SEB derives the latent heat flux lazily during the skin
    # temperature solve (see the lazy `surface_humidity_flux` method below).
    out.canopy_evaporation_conductance[i, j, 1] = canopy_evaporation_conductance(evapotranspiration, f_can, rₐ)
    out.transpiration_conductance[i, j, 1] = transpiration_conductance(evapotranspiration, rₐ, gw_can)

    # Compute and store the partitioned ET fluxes from the current skin temperature. The ground
    # evaporation is driven by the ground temperature; the canopy evaporation and transpiration are
    # skin-driven. Re-running this kernel after the SEB solve (see `LandModel`) refreshes the
    # skin-driven fluxes so they are consistent with the converged skin temperature.
    out.evaporation_ground[i, j, 1] = compute_evaporation_ground(evapotranspiration, Δqg, β, rₐ, rₑ)
    out.transpiration[i, j, 1] = compute_transpiration(evapotranspiration, Δqs, rₐ, gw_can)
    out.evaporation_canopy[i, j, 1] = compute_evaporation_canopy(evapotranspiration, Δqs, f_can, rₐ)
    return out
end

"""
    $TYPEDSIGNATURES

Compute the aerodynamic resistance between the ground and canopy as a function of LAI and SAI.
"""
@inline function aerodynamic_resistance(
        i, j, grid, fields,
        atmos::AbstractAtmosphere,
        evapotranspiration::PALADYNCanopyEvapotranspiration{NF},
        ::AbstractVegetation # included just to make explicit the dependence on vegetation fields
    ) where {NF}
    @inbounds let LAI = max(fields.leaf_area_index[i, j], zero(NF)),
            SAI = max(fields.stem_area_index[i, j], zero(NF)),
            Vₐ = windspeed(i, j, grid, fields, atmos),
            C = evapotranspiration.C_can  # drag coefficient for the canopy
        rₙ = (1 - exp(-LAI - SAI)) / (C * Vₐ)
        return rₙ
    end
end
