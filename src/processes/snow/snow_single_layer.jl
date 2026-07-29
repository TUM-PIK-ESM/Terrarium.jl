"""
    $TYPEDEF

Simple single-layer snow scheme (loosely based on the Utah Energy Balance model,
[tarbotonSpatiallyDistributedEnergy1994](@cite)). The snowpack is represented as a single lumped layer
with a bulk density `ρ_snow` supplied by a snow-density scheme (constant by default), from which the thermal
properties follow. The prognostic state is the depth-integrated (column) internal energy `snow_energy`
`Ū_snow` [J/m²] and the `snow_water_equivalent` `W_snow` [m]; snow depth, cover fraction, and thermal conductivity
are diagnosed from these and the bulk density.

Properties:
$FIELDS

# References
* [tarbotonSpatiallyDistributedEnergy1994](@cite) Tarboton, Chowdhury and Jackson (1994)
"""
@parameterized struct SingleLayerSnow{NF, Cover, Density, Conductivity, Hydraulics, Closure} <: AbstractSnow{NF}
    "Snow areal coverage parameterization"
    @component cover::Cover

    "Bulk snow density parameterization"
    @component density::Density

    "Snow thermal conductivity parameterization"
    @component thermal_conductivity::Conductivity

    "Snow hydraulic properties"
    @component hydraulic_properties::Hydraulics

    "Snow energy-temperature closure"
    @component closure::Closure
end

function SingleLayerSnow(
        ::Type{NF};
        cover = FractionalSnowCover(NF),
        density = ConstantSnowDensity(NF),
        thermal_conductivity = PowerLawSnowThermalConductivity(NF),
        hydraulic_properties = ConstantSnowHydraulics(NF),
        closure = SnowEnergyTemperatureClosure(NF)
    ) where {NF}
    # `NF` is not carried by any field, so it must be supplied explicitly to the type constructor
    return SingleLayerSnow{NF, typeof(cover), typeof(density), typeof(thermal_conductivity), typeof(hydraulic_properties), typeof(closure)}(
        cover, density, thermal_conductivity, hydraulic_properties, closure
    )
end

# Top-level interface

"""
    $TYPEDSIGNATURES

Snow layer depth `d_snow = W_snow·ρ_w/ρ_snow` [m], converting the water-equivalent depth `W_snow` [m] to the physical
snow depth using the water density `ρ_w` and the bulk snow density `ρ_snow`.
"""
@inline compute_snow_depth(::AbstractSnow, W_snow::NF, ρ_snow::NF, ρ_w::NF) where {NF} = max(W_snow, zero(NF)) * ρ_w / ρ_snow

"""
    $TYPEDSIGNATURES

Bulk snow density `ρ_snow` [kg/m³] of the snowpack, delegating to the process's density scheme.
"""
@inline snow_density(snow::SingleLayerSnow) = snow_density(snow.density)

"""
    $TYPEDSIGNATURES

Bulk snow thermal conductivity `κ_snow` [W/m/K], delegating to the process's thermal conductivity scheme
with the bulk density `ρ_snow`.
"""
@inline compute_thermal_conductivity(snow::SingleLayerSnow, constants::MaterialConstants, ρ_snow) =
    compute_thermal_conductivity(snow.thermal_conductivity, constants, ρ_snow)

"""$TYPEDSIGNATURES"""
@inline get_closure(snow::SingleLayerSnow) = snow.closure

"""
    $TYPEDSIGNATURES

Launch [`compute_snow_soil_boundary_fluxes!`](@ref) to diagnose the snow↔surface/soil coupling fluxes from the
snow state and the surface energy balance outputs: the blended soil-top heat flux (`soil_heat_flux`) and
the snow surface sublimation rate (`sublimation`, the snow-fraction bulk-aerodynamic vapor flux at the
converged skin temperature). Must run after the surface energy balance (which sets `ground_heat_flux` and
the skin temperature). No-op when there is no snowpack (`snow === nothing`).
"""
compute_snow_soil_boundary_fluxes!(state, grid, ::Nothing, constants::PhysicalConstants, atmos::AbstractAtmosphere) = nothing
function compute_snow_soil_boundary_fluxes!(state, grid, snow::SingleLayerSnow, constants::PhysicalConstants, atmos::AbstractAtmosphere)
    out = (; soil_heat_flux = state.soil_heat_flux, sublimation = state.sublimation)
    coupling_fields = (;
        ground_heat_flux = state.ground_heat_flux,
        ground_temperature = state.ground_temperature,
        skin_temperature = state.skin_temperature,
        snow_temperature = state.snow_temperature,
        snow_depth = state.snow_depth,
        snow_cover_fraction = state.snow_cover_fraction,
    )
    # atmosphere fields (air temperature/pressure/humidity, windspeed) are needed for the sublimation flux
    fields = merge(coupling_fields, get_fields(state, atmos))
    launch!(grid, XY, compute_snow_coupling_fluxes_kernel!, out, fields, snow, atmos, constants)
    return nothing
end

# Process methods

variables(snow::SingleLayerSnow) = (
    prognostic(:snow_energy, XY(); closure = get_closure(snow), units = u"J/m^2", desc = "Depth-integrated (column) internal energy of the snowpack relative to water at 0°C"),
    prognostic(:snow_water_equivalent, XY(); units = u"m", desc = "Snow water equivalent (ice + retained liquid)"),
    auxiliary(:snow_depth, XY(); units = u"m", desc = "Snow layer depth"),
    auxiliary(:snow_cover_fraction, XY(); bounds = UnitInterval, desc = "Sub-grid snow-covered area fraction"),
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

Accumulate the snow water-equivalent and depth-integrated energy tendencies at grid cell `i, j` from the
mass and energy balances (see [`compute_snow_water_tendency`](@ref) and [`compute_snow_energy_tendency`](@ref)).
"""
@propagate_inbounds function compute_snow_tendencies!(
        tendencies, i, j, grid, fields,
        snow::SingleLayerSnow,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants
    )
    tendencies.snow_water_equivalent[i, j, 1] += compute_snow_water_tendency(i, j, grid, fields, snow, atmos)
    tendencies.snow_energy[i, j, 1] += compute_snow_energy_tendency(i, j, grid, fields, snow, atmos, constants)
    return nothing
end
"""
    $TYPEDSIGNATURES

Compute the snow depth, cover fraction, and thermal conductivity at grid cell `i, j`.
"""
@propagate_inbounds function compute_snow_properties!(
        out, i, j, grid, fields,
        snow::SingleLayerSnow,
        constants::PhysicalConstants
    )
    W_snow = fields.snow_water_equivalent[i, j]
    ρ_w = constants.material.density_water
    ρ_snow = compute_snow_density(i, j, grid, fields, snow.density)
    out.snow_depth[i, j, 1] = compute_snow_depth(snow, W_snow, ρ_snow, ρ_w)
    out.snow_cover_fraction[i, j, 1] = compute_snow_cover_fraction(snow.cover, W_snow)
    return nothing
end

"""
    $TYPEDSIGNATURES

Blended soil-top heat flux [W/m²] at grid cell `i, j`: the snow-cover-fraction-weighted combination of
the snow→soil basal conductive flux `Q_base` and the bare-ground surface energy balance closure flux `G`
(`ground_heat_flux`), `f_snow·Q_base + (1 − f_snow)·G`. The bulk snow thermal conductivity entering
`Q_base` is recovered lazily from the density scheme rather than stored.
"""
@propagate_inbounds function compute_snow_soil_heat_flux(i, j, grid, fields, snow::SingleLayerSnow, constants::PhysicalConstants)
    G = fields.ground_heat_flux[i, j]
    f = fields.snow_cover_fraction[i, j]
    d_snow = fields.snow_depth[i, j]
    T_snow = fields.snow_temperature[i, j]
    T_soil = fields.ground_temperature[i, j]
    ρ_snow = compute_snow_density(i, j, grid, fields, snow.density)
    κ_snow = compute_thermal_conductivity(snow, constants.material, ρ_snow)
    Q_base = compute_snow_basal_heat_flux(κ_snow, T_soil, T_snow, d_snow)
    return f * Q_base + (one(f) - f) * G
end

"""
    $TYPEDSIGNATURES

Diagnose the snow↔surface/soil coupling fluxes at grid cell `i, j`, run after the surface energy balance:
the blended soil-top heat flux (see [`compute_snow_soil_heat_flux`](@ref)) and the snow surface
sublimation rate (see [`compute_snow_sublimation_flux`](@ref), the same snow-fraction vapor flux the
surface energy balance uses for the latent-flux partition). The snow surface and basal heat fluxes seen
by the energy tendency are the `ground_heat_flux` and `soil_heat_flux` fields themselves (aliased as the
snow's `surface_heat_flux`/`basal_heat_flux` inputs in the coupled model).
"""
@propagate_inbounds function compute_snow_soil_boundary_fluxes!(
        out, i, j, grid, fields,
        snow::SingleLayerSnow,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants
    )
    out.soil_heat_flux[i, j, 1] = compute_snow_soil_heat_flux(i, j, grid, fields, snow, constants)
    out.sublimation[i, j, 1] = compute_snow_sublimation_flux(i, j, grid, fields, snow, atmos, constants)
    return nothing
end

"""
    $TYPEDSIGNATURES

Snow meltwater outflow `M_r` [m/s SWE] at grid cell `i, j`: the Darcy-type drainage
(see [`compute_meltwater_outflow`](@ref)) evaluated at the diagnosed liquid water fraction.
"""
@propagate_inbounds function snow_meltwater_flux(i, j, grid, fields, snow::SingleLayerSnow)
    θ_liq = liquid_water_fraction(i, j, grid, fields, snow)
    M = compute_meltwater_outflow(snow.hydraulic_properties, θ_liq)
    return M
end

## Field-level accessors

@propagate_inbounds snow_water_equivalent(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_water_equivalent[i, j]

@propagate_inbounds snow_energy(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_energy[i, j]

@propagate_inbounds snow_depth(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_depth[i, j]

@propagate_inbounds snow_cover_fraction(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_cover_fraction[i, j]

@propagate_inbounds snow_temperature(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_temperature[i, j]

@propagate_inbounds liquid_water_fraction(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_liquid_fraction[i, j]

# Kernels

@kernel inbounds = true function compute_snow_coupling_fluxes_kernel!(out, grid, fields, snow::SingleLayerSnow, atmos, constants)
    i, j = @index(Global, NTuple)
    compute_snow_soil_boundary_fluxes!(out, i, j, grid, fields, snow, atmos, constants)
end

@kernel inbounds = true function compute_auxiliary_kernel!(out, grid, fields, snow::AbstractSnow, args...)
    i, j = @index(Global, NTuple)
    compute_snow_properties!(out, i, j, grid, fields, snow, args...)
end

@kernel inbounds = true function compute_tendencies_kernel!(tendencies, grid, fields, snow::SingleLayerSnow, args...)
    i, j = @index(Global, NTuple)
    compute_snow_tendencies!(tendencies, i, j, grid, fields, snow, args...)
end
