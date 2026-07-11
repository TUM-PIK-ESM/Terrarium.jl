"""
    $TYPEDEF

Inert (no-op) snow component representing the absence of a snowpack. Analogous to
[`NoCanopyInterception`](@ref); the default land-model configuration uses `NoSnow` so that non-snow
behavior is unchanged. All snow accessors return zero.
"""
struct NoSnow{NF} <: AbstractSnow{NF} end

NoSnow(::Type{NF}) where {NF} = NoSnow{NF}()

variables(::NoSnow) = ()

@inline compute_auxiliary!(state, grid, ::NoSnow, args...) = nothing

@inline compute_tendencies!(state, grid, ::NoSnow, args...) = nothing

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
@parameterized @kwdef struct SingleLayerSnow{NF, Density} <: AbstractSnow{NF}
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
    @param density::Density = ConstantSnowDensity(typeof(cover_reference))
end

SingleLayerSnow(::Type{NF}; density = ConstantSnowDensity(NF), kwargs...) where {NF} = SingleLayerSnow{NF, typeof(density)}(; density, kwargs...)

variables(snow::SingleLayerSnow) = (
    prognostic(:snow_energy, XY(); closure = get_closure(snow), units = u"J/m^2", desc = "Depth-integrated (column) internal energy of the snowpack relative to ice at 0°C"),
    prognostic(:snow_water_equivalent, XY(); units = u"m", desc = "Snow water equivalent (ice + retained liquid)"),
    auxiliary(:snow_depth, XY(); units = u"m", desc = "Snow layer depth"),
    auxiliary(:snow_cover_fraction, XY(); bounds = UnitInterval, desc = "Sub-grid snow-covered area fraction"),
    auxiliary(:snow_thermal_conductivity, XY(); units = u"W/m/K", desc = "Bulk snow thermal conductivity"),
    # Surface forcing consumed by the energy/mass tendencies. In coupling these are provided by the
    # surface energy balance (`surface_heat_flux`, `sublimation`) and the snow→soil conduction
    # (`basal_heat_flux`); the standalone `SnowModel` prescribes them directly.
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
    out.snow_depth[i, j, 1] = compute_snow_depth(snow, W, ρ_w)
    out.snow_cover_fraction[i, j, 1] = compute_snow_cover_fraction(snow, W)
    out.snow_thermal_conductivity[i, j, 1] = compute_snow_thermal_conductivity(snow, ρ_w)
    return nothing
end

# Kernels

@kernel inbounds = true function compute_auxiliary_kernel!(out, grid, fields, snow::AbstractSnow, args...)
    i, j = @index(Global, NTuple)
    compute_snow_properties!(out, i, j, grid, fields, snow, args...)
end
