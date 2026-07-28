"""
    $TYPEDEF

Energy–temperature closure for the single-layer snow scheme. Recovers the depth-averaged snow
temperature `T_snow` [°C] and liquid water fraction `θ_liq` from the depth-integrated internal energy
`E` [J/m²] using the medium-agnostic `FreeWater` enthalpy relations, treating the bulk snow as a water
substance of volumetric mass `ρ_s`: the volumetric energy `U_v = E/d_s` [J/m³], the volumetric latent
heat `Lθ = ρ_s·L_{sl}`, and the volumetric heat capacity `C = ρ_s·(θ_liq·c_{p,l} + (1−θ_liq)·c_{p,i})`. The
volumetric energy `U_v` is formed inline from the diagnosed snow depth `d_s` (no auxiliary field is
stored for `U_v`). This is the same enthalpy map used for soil with the pore saturation set to one,
i.e. the snow layer is "fully saturated" with water substance in the lumped sense.
"""
struct SnowEnergyTemperatureClosure <: AbstractEnergyClosure end

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

"""
    $TYPEDSIGNATURES

Inverse closure: compute the prognostic `snow_energy` from a prescribed `snow_temperature` (used for
initialization). `snow_water_equivalent` must be initialized first, since the snow depth `d_s` enters
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
    E = fields.snow_energy[i, j] # assumed given
    W = fields.snow_water_equivalent[i, j]
    L_sl = constants.thermodynamics.latent_heat_fusion
    c_pi = constants.thermodynamics.specific_heat_capacity_ice
    c_pw = constants.thermodynamics.specific_heat_capacity_liquid_water
    ρ_w = constants.material.density_water
    ρ_s = compute_snow_density(i, j, grid, fields, snow.density)
    d_s = compute_snow_depth(snow, W, ρ_s, ρ_w)
    # Volumetric latent heat of fusion and volumetric energy (bulk snow treated as water substance)
    Lθ = ρ_s * L_sl
    U_v = volumetric_snow_energy(E, d_s)
    liq = liquid_water_fraction(FreeWater(), U_v, Lθ, one(U_v))
    out.snow_liquid_fraction[i, j, 1] = liq
    C = ρ_s * (liq * c_pw + (one(liq) - liq) * c_pi)
    # Snow temperature cannot exceed 0°C, so clip the free-water temperature at zero. The energy above
    # the fully-melted (0°C, all-liquid) reference, i.e. the positive part `U_v > 0`, is not stored; it
    # is derived on demand where needed to determine snow melt.
    T = energy_to_temperature(FreeWater(), U_v, Lθ, C)
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
    T = fields.snow_temperature[i, j] # assumed given
    W = fields.snow_water_equivalent[i, j]
    L_sl = constants.thermodynamics.latent_heat_fusion
    c_pi = constants.thermodynamics.specific_heat_capacity_ice
    c_pw = constants.thermodynamics.specific_heat_capacity_liquid_water
    ρ_w = constants.material.density_water
    ρ_s = compute_snow_density(i, j, grid, fields, snow.density)
    d_s = compute_snow_depth(snow, W, ρ_s, ρ_w)
    Lθ = ρ_s * L_sl
    # N.B. For the free-water characteristic the liquid fraction is indeterminate at T = 0; assume
    # frozen for T < 0 and thawed otherwise. This mapping is for initialization only and must **not**
    # be used in the calculation of tendencies.
    liq = ifelse(T >= zero(T), one(T), zero(T))
    out.snow_liquid_fraction[i, j, 1] = liq
    C = ρ_s * (liq * c_pw + (one(liq) - liq) * c_pi)
    U_v = T * C - Lθ * (one(liq) - liq)
    out.snow_energy[i, j, 1] = U_v * d_s
    return nothing
end

"""
    $TYPEDSIGNATURES

Snow→soil basal conductive heat flux `Q_base` [W/m²], positive upward (soil → snow). The snowpack is a
strong insulator, so only its resistance is retained (snow-resistance-only closure):
`Q_base = 2·κ_snow·(T_soil − T_snow)/d_s`. The depth is regularized with a machine-`eps` offset so the
flux stays finite as the snowpack vanishes.
"""
@inline compute_snow_basal_heat_flux(κ_snow::NF, T_soil::NF, T_snow::NF, d_s::NF) where {NF} =
    2 * κ_snow * (T_soil - T_snow) / (d_s + eps(NF))

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

# Helper functions

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
