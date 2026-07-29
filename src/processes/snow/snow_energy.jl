"""
    $TYPEDSIGNATURES

Advected heat flux [W/m²] carried into the snowpack by precipitation, relative to liquid water at 0 °C
(the `U = 0` reference of the `FreeWater` enthalpy closure). Fresh snow `P_s` arrives as ice, which sits
`L_f` below the liquid reference, plus sensible heat for `T_air < 0`; rain-on-snow `R_on_snow` arrives as
liquid carrying only its sensible heat for `T_air > 0`. The latent heat released when rain refreezes in a
cold pack is captured implicitly by the enthalpy closure, so it is *not* added here (adding `L_f` would
double-count relative to the liquid-water reference).
"""
@inline function compute_snow_precip_heat_flux(::AbstractSnow, constants::PhysicalConstants, P_s::NF, R_on_snow::NF, T_air::NF) where {NF}
    ρ_w = constants.material.density_water
    L_sl = constants.thermodynamics.latent_heat_fusion
    cp_i = constants.thermodynamics.specific_heat_capacity_ice
    cp_w = constants.thermodynamics.specific_heat_capacity_liquid_water
    Q_snow = ρ_w * P_s * (cp_i * min(T_air, zero(NF)) - L_sl)
    Q_rain = ρ_w * R_on_snow * cp_w * max(T_air, zero(NF))
    Q_prcp = Q_snow + Q_rain
    return Q_prcp
end

"""
    $TYPEDSIGNATURES

Bulk volumetric heat capacity `C_snow` [J/m³/K] of the snowpack treated as an ice–liquid–air mixture,
given the bulk snow density `ρ_snow` [kg/m³] and the liquid water fraction `liq` ∈ [0,1] of the water
substance. The water-substance mass per unit snow volume is `ρ_snow`, of which `(1 − liq)` is ice and
`liq` is liquid; the remaining void space is dry air. The corresponding constituent volume fractions are
```math
\\begin{aligned}
θ_ice &= ρ_snow·(1 − liq)/ρ_i,\\
θ_liq &= ρ_snow·liq/ρ_w,\\
θ_air &= 1 − θ_ice − θ_liq
\\end{aligned}
```
and the heat capacity is the volume-weighted sum over the ice and liquid constituents,
`C_snow = cp_i·ρ_i·θ_ice + cp_w·ρ_w·θ_liq`. The air's own sensible-heat storage (≈0.1%) is neglected; note
the air is still reflected in the bulk density `ρ_snow < ρ_ice`.
"""
@propagate_inbounds function compute_snow_volumetric_heat_capacity(
        ::AbstractSnow{NF},
        constants::PhysicalConstants,
        ρ_snow::NF,
        liq::NF
    ) where {NF}
    cp_i = constants.thermodynamics.specific_heat_capacity_ice
    cp_w = constants.thermodynamics.specific_heat_capacity_liquid_water
    ρ_w = constants.material.density_water
    ρ_i = constants.material.density_ice
    # ice and liquid volume fractions per unit snow volume; the remaining void is air, whose
    # sensible-heat storage (≈0.1%) is neglected
    θ_ice = ρ_snow * (one(NF) - liq) / ρ_i
    θ_liq = ρ_snow * liq / ρ_w
    C_snow = cp_i * ρ_i * θ_ice + cp_w * ρ_w * θ_liq
    return C_snow
end

"""
    $TYPEDSIGNATURES

Depth-integrated snow energy tendency [W/m²] at grid cell `i, j` (all fluxes positive upward):
```
dŪ_snow/dt = Q_base − Q_top + Q_precip + Q_subl
```
where `Q_top`/`Q_base` are the surface/basal heat fluxes, `Q_precip` the advected precipitation heat
(see [`compute_snow_precip_heat_flux`](@ref)), and `Q_subl` an advective correction for sublimation.

The sublimation correction `Q_subl = ρ_w·L_f·E_subl` is required because the latent heat flux
carries the full sublimation enthalpy `ρ_w·L_sg·E_subl`, whereas the mass leaving the snowpack departs as ice,
whose specific enthalpy relative to the liquid-water reference is `−L_f`. Adding back `ρ_w·L_f·E_subl`
leaves the snowpack with a net loss of `ρ_w·(L_s − L_f)·E_subl = ρ_w·L_v·E_subl`, the vaporization enthalpy
carried by the departing vapor.

Note that no explicit *meltwater* energy term appears because meltwater drains as liquid water at 0 °C,
which is the zero-enthalpy reference (`U = 0`) of the `FreeWater` closure, so it carries no enthalpy
out of the snowpack.
"""
@propagate_inbounds function compute_snow_energy_tendency(
        i, j, grid, fields,
        snow::SingleLayerSnow,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants
    )
    ρ_w = constants.material.density_water
    L_sl = constants.thermodynamics.latent_heat_fusion
    Q_top = fields.surface_heat_flux[i, j]   # positive upward: energy leaving the snow top
    Q_base = fields.basal_heat_flux[i, j]    # positive upward: energy entering the snow base from soil
    E_subl = fields.sublimation[i, j]        # positive upward: snow mass loss to water vapor
    f_snow = snow_cover_fraction(i, j, grid, fields, snow)
    P_s = snowfall(i, j, grid, fields, atmos)
    R = rainfall(i, j, grid, fields, atmos)
    T_air = air_temperature(i, j, grid, fields, atmos)
    R_snow = f_snow * R
    Q_prcp = compute_snow_precip_heat_flux(snow, constants, P_s, R_snow, T_air)
    Q_subl = ρ_w * L_sl * E_subl # correction to account for enthalpy of lost ice
    dUdt = Q_base - Q_top + Q_prcp + Q_subl
    return dUdt
end

# Helper functions

"""
    $TYPEDSIGNATURES

Volumetric snow internal energy `U_snow = Ū_snow/d_snow` [J/m³] from the depth-integrated energy `Ū_snow`
[J/m²] and the snow depth `d_snow` [m]. The denominator is regularized with a machine-`eps` offset so the
result stays finite (and differentiable) as `W_snow → 0`: with `Ū_snow → 0` and `d_snow → 0` together,
`U_snow → 0`. The thin-snow indeterminacy is additionally masked downstream by `f_snow → 0`. A finite
offset is used rather than [`safediv`](@ref) because `safediv` returns `Inf` at exactly `d_snow = 0` (even
for `Ū_snow = 0`), which would produce `NaN` when multiplied by `f_snow = 0` in the surface energy
balance blend.
"""
@inline compute_volumetric_snow_energy(Ū_snow::NF, d_snow::NF) where {NF} = Ū_snow / (d_snow + eps(NF))

"""
    $TYPEDSIGNATURES

Snow→soil basal conductive heat flux `Q_base` [W/m²], positive upward (soil → snow). The snowpack is a
strong insulator, so only its resistance is retained (snow-resistance-only closure):
`Q_base = 2·κ_snow·(T_soil − T_snow)/d_snow`. The depth is regularized with a machine-`eps` offset so the
flux stays finite as the snowpack vanishes.
"""
@inline compute_snow_basal_heat_flux(κ_snow::NF, T_soil::NF, T_snow::NF, d_snow::NF) where {NF} =
    2 * κ_snow * (T_soil - T_snow) / (d_snow + eps(NF))


# Snow energy closure

"""
    $TYPEDEF

Energy–temperature closure for snow volumes. For [`SingleLayerSnow`](@ref), the depth-averaged snow
temperature `T_snow` (°C) and liquid water fraction `θ_liq` are recovered from the depth-integrated internal energy
`Ū_snow` (J/m²) using the medium-agnostic `FreeWater` enthalpy relations, treating the bulk snowpack as an ice-water-air
mixture.
"""
struct SnowEnergyTemperatureClosure{NF} <: AbstractEnergyClosure end

SnowEnergyTemperatureClosure(::Type{NF}) where {NF} = SnowEnergyTemperatureClosure{NF}()

variables(::SnowEnergyTemperatureClosure) = (
    auxiliary(:snow_temperature, XY(), units = u"°C", desc = "Depth-averaged snow temperature in °C (≤ 0)"),
    auxiliary(:snow_liquid_fraction, XY(), bounds = UnitInterval, desc = "Liquid (unfrozen) fraction of the snow water substance"),
)

# Process-level closure entry points (dispatch on the snow process, mirroring the soil interface).
@inline closure!(state, grid, snow::SingleLayerSnow, constants::PhysicalConstants) =
    closure!(state, grid, get_closure(snow), snow, constants)

@inline invclosure!(state, grid, snow::SingleLayerSnow, constants::PhysicalConstants) =
    invclosure!(state, grid, get_closure(snow), snow, constants)

"""
    $TYPEDSIGNATURES

Initialize the snow internal energy by evaluating the inverse closure (temperature → energy). Assumes
`snow_temperature` and `snow_water_equivalent` have already been initialized.
"""
function initialize!(state, grid, snow::SingleLayerSnow, constants::PhysicalConstants, args...)
    invclosure!(state, grid, get_closure(snow), snow, constants)
    return nothing
end

"""
    $TYPEDSIGNATURES

Forward closure: recover `snow_temperature` and `snow_liquid_fraction` from the prognostic `snow_energy`.
"""
function closure!(
        state, grid,
        closure::SnowEnergyTemperatureClosure,
        snow::SingleLayerSnow,
        constants::PhysicalConstants,
        args...
    )
    kernel_args = (closure, snow, constants)
    out = closure_fields(state, snow)
    fields = get_fields(state, kernel_args...; except = out)
    launch!(grid, XY, snow_energy_to_temperature_kernel!, out, fields, kernel_args...)
    return nothing
end

# Helper functions

"""
    $TYPEDSIGNATURES

Inverse closure: compute the prognostic `snow_energy` from a prescribed `snow_temperature` (used for
initialization). `snow_water_equivalent` must be initialized first, since the snow depth `d_snow` enters
the depth-integrated energy.
"""
function invclosure!(
        state, grid,
        closure::SnowEnergyTemperatureClosure,
        snow::SingleLayerSnow,
        constants::PhysicalConstants,
        args...
    )
    kernel_args = (closure, snow, constants)
    # snow_energy is the prognostic variable, so collect the output fields manually
    out = (snow_energy = state.snow_energy, snow_liquid_fraction = state.snow_liquid_fraction)
    fields = get_fields(state, kernel_args...; except = out)
    launch!(grid, XY, snow_temperature_to_energy_kernel!, out, fields, kernel_args...)
    return nothing
end

# Kernel functions

"""
    $TYPEDSIGNATURES

Recover the snow temperature and liquid water fraction from the depth-integrated energy at grid cell `i, j`.
"""
@propagate_inbounds function energy_to_temperature!(
        out, i, j, grid, fields,
        ::SnowEnergyTemperatureClosure,
        snow::SingleLayerSnow,
        constants::PhysicalConstants
    )
    Ū_snow = fields.snow_energy[i, j] # assumed given (prognostic)
    W_snow = fields.snow_water_equivalent[i, j] # assumed given (prognostic)
    L_sl = constants.thermodynamics.latent_heat_fusion
    ρ_w = constants.material.density_water
    ρ_snow = compute_snow_density(i, j, grid, fields, snow.density)
    d_snow = compute_snow_depth(snow, W_snow, ρ_snow, ρ_w)
    # Volumetric latent heat of fusion
    Lθ = ρ_snow * L_sl
    U_snow = compute_volumetric_snow_energy(Ū_snow, d_snow)
    liq = liquid_water_fraction(FreeWater(), U_snow, Lθ)
    out.snow_liquid_fraction[i, j, 1] = liq
    C_snow = compute_snow_volumetric_heat_capacity(snow, constants, ρ_snow, liq)
    # Snow temperature cannot exceed 0°C, so clip the free-water temperature at zero. The energy above
    # the fully-melted (0°C, all-liquid) reference, i.e. the positive part `U_snow > 0`, is not stored; it
    # is derived on demand where needed to determine snow melt.
    T = energy_to_temperature(FreeWater(), U_snow, Lθ, C_snow)
    out.snow_temperature[i, j, 1] = min(T, zero(T))
    return nothing
end

"""
    $TYPEDSIGNATURES

Compute the depth-integrated snow energy from a prescribed temperature at grid cell `i, j`.
"""
@propagate_inbounds function temperature_to_energy!(
        out, i, j, grid, fields,
        ::SnowEnergyTemperatureClosure,
        snow::SingleLayerSnow,
        constants::PhysicalConstants
    )
    T = fields.snow_temperature[i, j] # assumed given (initialization)
    W_snow = fields.snow_water_equivalent[i, j] # assumed given (prognostic)
    L_sl = constants.thermodynamics.latent_heat_fusion
    ρ_w = constants.material.density_water
    ρ_snow = compute_snow_density(i, j, grid, fields, snow.density)
    d_snow = compute_snow_depth(snow, W_snow, ρ_snow, ρ_w)
    Lθ = ρ_snow * L_sl
    # N.B. For the free-water characteristic the liquid fraction is indeterminate at T = 0; assume
    # frozen for T < 0 and thawed otherwise. This mapping is for initialization only and must **not**
    # be used in the calculation of tendencies.
    liq = ifelse(T >= zero(T), one(T), zero(T))
    out.snow_liquid_fraction[i, j, 1] = liq
    C_snow = compute_snow_volumetric_heat_capacity(snow, constants, ρ_snow, liq)
    U_snow = T * C_snow - Lθ * (one(liq) - liq)
    out.snow_energy[i, j, 1] = U_snow * d_snow # integrate by d_snow
    return nothing
end

# Kernels
#
# These wrap the shared `energy_to_temperature!` / `temperature_to_energy!` kernel functions (dispatched
# on the snow closure type) in 2D (XY) kernels. They use snow-specific names because KernelAbstractions
# `@kernel` does not support multiple methods of one kernel — the soil energy closures define 3D (XYZ)
# kernels under the base names.

@kernel inbounds = true function snow_energy_to_temperature_kernel!(out, grid, fields, closure::SnowEnergyTemperatureClosure, args...)
    i, j = @index(Global, NTuple)
    energy_to_temperature!(out, i, j, grid, fields, closure, args...)
end

@kernel inbounds = true function snow_temperature_to_energy_kernel!(out, grid, fields, closure::SnowEnergyTemperatureClosure, args...)
    i, j = @index(Global, NTuple)
    temperature_to_energy!(out, i, j, grid, fields, closure, args...)
end
