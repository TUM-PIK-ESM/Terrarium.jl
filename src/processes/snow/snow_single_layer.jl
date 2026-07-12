"""
    $TYPEDEF

Simple single-layer snow scheme (loosely based on the Utah Energy Balance model,
[tarbotonSpatiallyDistributedEnergy1994](@cite)). The snowpack is represented as a single lumped layer
with a bulk density `ρ_s` supplied by a snow-density scheme (constant by default), from which the thermal
properties follow. The prognostic state is the depth-integrated (column) internal energy `snow_energy`
`E` [J/m²] and the `snow_water_equivalent` `W` [m]; snow depth, cover fraction, and thermal conductivity
are diagnosed from these and the bulk density.

Properties:
$FIELDS

# References
* [tarbotonSpatiallyDistributedEnergy1994](@cite) Tarboton, Chowdhury and Jackson (1994)
"""
@parameterized @kwdef struct SingleLayerSnow{NF, Density, Closure} <: AbstractSnow{NF}
    "Reference snow water equivalent `W_ref` for the fractional snow-cover diagnostic `f_snow = W/(W + W_ref)`"
    @param cover_reference::NF = 0.01 (units = u"m", bounds = Positive)

    "Coefficient `a` in the thermal conductivity power law `κ = a·(ρ_s/ρ_w)^b` ([yen1981review](@cite))"
    @param conductivity_coefficient::NF = 2.22362 (units = u"W/m/K", bounds = Positive)

    "Exponent `b` in the thermal conductivity power law `κ = a·(ρ_s/ρ_w)^b` ([yen1981review](@cite))"
    @param conductivity_exponent::NF = 1.885 (bounds = Positive,)

    "Capillary retention `L_c`: liquid fraction held against gravity before meltwater drains"
    @param capillary_retention::NF = 0.05 (bounds = UnitInterval)

    "Saturated hydraulic conductivity `K_sat` of the snowpack, setting the meltwater outflow rate"
    @param saturated_conductivity::NF = 1.0e-4 (units = u"m/s", bounds = Positive)

    "Bulk snow density scheme"
    @component density::Density = ConstantSnowDensity(typeof(cover_reference))

    "Snow energy-temperature closure"
    @component closure::Closure = SnowEnergyTemperatureClosure()
end

function SingleLayerSnow(
        ::Type{NF};
        density = ConstantSnowDensity(NF),
        closure = SnowEnergyTemperatureClosure(),
        kwargs...
    ) where {NF}
    return SingleLayerSnow{NF, typeof(density), typeof(closure)}(; density, closure, kwargs...)
end

"""
    $TYPEDSIGNATURES

Snow layer depth `d_s = W·ρ_w/ρ_s` [m], converting the water-equivalent depth `W` [m] to the physical
snow depth using the water density `ρ_w` and the bulk snow density `ρ_s`.
"""
@inline compute_snow_depth(::AbstractSnow, W::NF, ρ_s::NF, ρ_w::NF) where {NF} = max(W, zero(NF)) * ρ_w / ρ_s

"""
    $TYPEDSIGNATURES

Sub-grid snow-covered area fraction `f_snow = W/(W + W_ref)` ∈ [0,1). Smooth and differentiable, with
`f_snow → 0` as `W → 0` and `f_snow → 1` as `W → ∞`.
"""
@inline function compute_snow_cover_fraction(snow::SingleLayerSnow, W::NF) where {NF}
    Wp = max(W, zero(NF))
    return Wp / (Wp + snow.cover_reference)
end

"""
    $TYPEDSIGNATURES

Bulk snow density `ρ_s` [kg/m³] of the snowpack, delegating to the process's density scheme.
"""
@inline snow_density(snow::SingleLayerSnow) = snow_density(snow.density)

"""
    $TYPEDSIGNATURES

Bulk snow thermal conductivity via the density power law `κ_snow = a·(ρ_s/ρ_w)^b`
([yen1981review](@cite)), constant for a fixed bulk density `ρ_s`.
"""
@inline compute_snow_thermal_conductivity(snow::SingleLayerSnow, ρ_w::NF) where {NF} =
    snow.conductivity_coefficient * (snow_density(snow) / ρ_w)^snow.conductivity_exponent

"""
    $TYPEDSIGNATURES

Volumetric snow internal energy `U_v = E/d_s` [J/m³] from the depth-integrated energy `E` [J/m²] and the
snow depth `d_s` [m]. The denominator is regularized with a machine-`eps` offset so the result stays
finite (and differentiable) as `W → 0`: with `E → 0` and `d_s → 0` together, `U_v → 0`. The thin-snow
indeterminacy is additionally masked downstream by `f_snow → 0`. A finite offset is used rather than
[`safediv`](@ref) because `safediv` returns `Inf` at exactly `d_s = 0` (even for `E = 0`), which would
produce `NaN` when multiplied by `f_snow = 0` in the surface energy balance blend.
"""
@inline volumetric_snow_energy(E::NF, d_s::NF) where {NF} = E / (d_s + eps(NF))

# Field-level accessors

@propagate_inbounds snow_water_equivalent(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_water_equivalent[i, j]

@propagate_inbounds snow_energy(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_energy[i, j]

@propagate_inbounds snow_depth(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_depth[i, j]

@propagate_inbounds snow_cover_fraction(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_cover_fraction[i, j]

@propagate_inbounds snow_temperature(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_temperature[i, j]

# Process methods

variables(snow::SingleLayerSnow) = (
    prognostic(:snow_energy, XY(); closure = get_closure(snow), units = u"J/m^2", desc = "Depth-integrated (column) internal energy of the snowpack relative to ice at 0°C"),
    prognostic(:snow_water_equivalent, XY(); units = u"m", desc = "Snow water equivalent (ice + retained liquid)"),
    auxiliary(:snow_depth, XY(); units = u"m", desc = "Snow layer depth"),
    auxiliary(:snow_cover_fraction, XY(); bounds = UnitInterval, desc = "Sub-grid snow-covered area fraction"),
    # Boundary heat fluxes consumed by the energy tendency. The standalone `SnowModel` prescribes them
    # directly; in the coupled `LandModel` they are aliased to the surface energy balance closure flux
    # (`surface_heat_flux ← ground_heat_flux`) and the blended soil-top flux (`basal_heat_flux ← soil_heat_flux`),
    # so no additional state fields are allocated (see the `LandModel` `StateVariables` constructor).
    input(:surface_heat_flux, XY(); units = u"W/m^2", desc = "Net heat flux at the snow surface (positive upward)"),
    input(:basal_heat_flux, XY(); units = u"W/m^2", desc = "Conductive heat flux at the snow base (positive upward, soil → snow)"),
    input(:sublimation, XY(); units = u"m/s", desc = "Sublimation/evaporation rate from the snow surface (SWE)"),
)

"""
    $TYPEDSIGNATURES

Diagnose the snow geometric and thermal properties (`snow_depth`, `snow_cover_fraction`,
`snow_thermal_conductivity`) from the current snow water equivalent and bulk density. The enthalpy
closure (`snow_temperature`, `snow_liquid_fraction`) is added in a later phase.
"""
function compute_auxiliary!(
        state, grid,
        snow::SingleLayerSnow,
        constants::PhysicalConstants,
        args...
    )
    out = auxiliary_fields(state, snow)
    fields = get_fields(state, snow; except = out)
    launch!(grid, XY, compute_auxiliary_kernel!, out, fields, snow, constants)
    return nothing
end

"""
    $TYPEDSIGNATURES

Compute the snow water-equivalent and depth-integrated energy tendencies (see [`compute_snow_tendencies!`](@ref)),
driven by the atmospheric precipitation inputs and the prescribed surface/basal heat fluxes and
sublimation.
"""
function compute_tendencies!(
        state, grid,
        snow::SingleLayerSnow,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        args...
    )
    tendencies = tendency_fields(state, snow)
    # pass the closure so `snow_liquid_fraction` (needed for meltwater outflow) is collected
    fields = get_fields(state, get_closure(snow), snow, atmos; except = tendencies)
    launch!(grid, XY, compute_tendencies_kernel!, tendencies, fields, snow, atmos, constants)
    return nothing
end

# Kernel functions

"""
    $TYPEDSIGNATURES

Compute the snow depth, cover fraction, and thermal conductivity at grid cell `i, j`.
"""
@propagate_inbounds function compute_snow_properties!(
        out, i, j, grid, fields,
        snow::SingleLayerSnow,
        constants::PhysicalConstants
    )
    W = fields.snow_water_equivalent[i, j]
    ρ_w = constants.material.density_water
    ρ_s = compute_snow_density(i, j, grid, fields, snow.density)
    out.snow_depth[i, j, 1] = compute_snow_depth(snow, W, ρ_s, ρ_w)
    out.snow_cover_fraction[i, j, 1] = compute_snow_cover_fraction(snow, W)
    return nothing
end

"""
    $TYPEDSIGNATURES

Diagnose the soil–snow interface heat flux at grid cell `i, j`, run after the surface energy balance.
The snow→soil basal conductive flux `Q_base` is computed lazily (the bulk thermal conductivity is
recovered from the snow-density scheme rather than stored), then blended into the soil-top boundary
flux by snow-covered area fraction: `soil_heat_flux = f_snow·Q_base + (1 − f_snow)·G`, where `G` is the
surface energy balance closure flux (`ground_heat_flux`). The snow surface and basal fluxes seen by the
energy tendency are the `ground_heat_flux` and `soil_heat_flux` fields themselves (aliased as the snow's
`surface_heat_flux`/`basal_heat_flux` inputs in the coupled model), so only `soil_heat_flux` is written here.
"""
@propagate_inbounds function compute_soil_snow_fluxes!(
        out, i, j, grid, fields,
        snow::SingleLayerSnow,
        constants::PhysicalConstants
    )
    ρ_w = constants.material.density_water
    κ_snow = compute_snow_thermal_conductivity(snow, ρ_w)
    G = fields.ground_heat_flux[i, j]
    f = fields.snow_cover_fraction[i, j]
    d_s = fields.snow_depth[i, j]
    T_snow = fields.snow_temperature[i, j]
    T_soil = fields.ground_temperature[i, j]
    Q_base = compute_snow_basal_heat_flux(κ_snow, T_soil, T_snow, d_s)
    # soil sees the fraction-weighted blend of the snow base (Q_base) and the bare ground (G)
    out.soil_heat_flux[i, j, 1] = f * Q_base + (one(f) - f) * G
    return nothing
end

# Kernels

@kernel inbounds = true function compute_auxiliary_kernel!(out, grid, fields, snow::AbstractSnow, args...)
    i, j = @index(Global, NTuple)
    compute_snow_properties!(out, i, j, grid, fields, snow, args...)
end

@kernel inbounds = true function compute_soil_snow_fluxes_kernel!(out, grid, fields, snow::SingleLayerSnow, constants)
    i, j = @index(Global, NTuple)
    compute_soil_snow_fluxes!(out, i, j, grid, fields, snow, constants)
end

# Top-level interface

"""
    $TYPEDSIGNATURES

Launch [`compute_soil_snow_fluxes!`](@ref) to diagnose the blended soil-top heat flux (`soil_heat_flux`)
from the snow state and the surface energy balance closure flux. Must run after the surface energy
balance (which sets `ground_heat_flux`). No-op when there is no snowpack (`snow === nothing`).
"""
compute_soil_snow_fluxes!(state, grid, ::Nothing, constants::PhysicalConstants) = nothing

function compute_soil_snow_fluxes!(state, grid, snow::SingleLayerSnow, constants::PhysicalConstants)
    out = (; soil_heat_flux = state.soil_heat_flux)
    fields = (;
        ground_heat_flux = state.ground_heat_flux,
        ground_temperature = state.ground_temperature,
        snow_temperature = state.snow_temperature,
        snow_depth = state.snow_depth,
        snow_cover_fraction = state.snow_cover_fraction,
    )
    launch!(grid, XY, compute_soil_snow_fluxes_kernel!, out, fields, snow, constants)
    return nothing
end
