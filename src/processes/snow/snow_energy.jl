"""
    $TYPEDEF

Energy–temperature closure for the single-layer snow scheme. Recovers the depth-averaged snow
temperature `T_snow` [°C] and liquid water fraction `θ_liq` from the depth-integrated internal energy
`E` [J/m²] using the medium-agnostic `FreeWater` enthalpy relations, treating the bulk snow as a water
substance of volumetric mass `ρ_s`: the volumetric energy `U_v = E/d_s` [J/m³], the volumetric latent
heat `Lθ = ρ_s·L_f`, and the volumetric heat capacity `C = ρ_s·(θ_liq·c_w + (1−θ_liq)·c_i)`. The
volumetric energy `U_v` is formed inline from the diagnosed snow depth `d_s` (no auxiliary field is
stored for `U_v`). This is the same enthalpy map used for soil with the pore saturation set to one,
i.e. the snow layer is "fully saturated" with water substance in the lumped sense.
"""
struct SnowEnergyTemperatureClosure <: AbstractEnergyClosure end

variables(::SnowEnergyTemperatureClosure) = (
    auxiliary(:snow_temperature, XY(), units = u"°C", desc = "Depth-averaged snow temperature in °C (≤ 0)"),
    auxiliary(:snow_liquid_fraction, XY(), bounds = UnitInterval, desc = "Liquid (unfrozen) fraction of the snow water substance"),
)

@inline get_closure(::SingleLayerSnow) = SnowEnergyTemperatureClosure()

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
    launch!(grid, XY, energy_to_temperature_kernel!, out, fields, kernel_args...)
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
    launch!(grid, XY, temperature_to_energy_kernel!, out, fields, kernel_args...)
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
    L_f = constants.thermodynamics.latent_heat_fusion
    c_i = constants.thermodynamics.specific_heat_capacity_ice
    c_w = constants.thermodynamics.specific_heat_capacity_liquid_water
    ρ_w = constants.material.density_water
    ρ_s = snow_density(i, j, grid, fields, snow.density)
    d_s = compute_snow_depth(snow, W, ρ_w)
    # Volumetric latent heat of fusion and volumetric energy (bulk snow treated as water substance)
    Lθ = ρ_s * L_f
    U_v = volumetric_snow_energy(E, d_s)
    liq = liquid_water_fraction(FreeWater(), U_v, Lθ, one(U_v))
    out.snow_liquid_fraction[i, j, 1] = liq
    C = ρ_s * (liq * c_w + (one(liq) - liq) * c_i)
    # Snow temperature cannot exceed 0°C, so clip the free-water temperature at zero. The energy above
    # the fully-melted (0°C, all-liquid) reference, i.e. the positive part `U_v > 0`, is not stored; it
    # is derived on demand where needed to drive melt and sublimation in the tendencies (later phase).
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
    L_f = constants.thermodynamics.latent_heat_fusion
    c_i = constants.thermodynamics.specific_heat_capacity_ice
    c_w = constants.thermodynamics.specific_heat_capacity_liquid_water
    ρ_w = constants.material.density_water
    ρ_s = snow_density(i, j, grid, fields, snow.density)
    d_s = compute_snow_depth(snow, W, ρ_w)
    Lθ = ρ_s * L_f
    # N.B. For the free-water characteristic the liquid fraction is indeterminate at T = 0; assume
    # frozen for T < 0 and thawed otherwise. This mapping is for initialization only and must **not**
    # be used in the calculation of tendencies.
    liq = ifelse(T >= zero(T), one(T), zero(T))
    out.snow_liquid_fraction[i, j, 1] = liq
    C = ρ_s * (liq * c_w + (one(liq) - liq) * c_i)
    U_v = T * C - Lθ * (one(liq) - liq)
    out.snow_energy[i, j, 1] = U_v * d_s
    return nothing
end

# Kernels
#
# These dispatch on the snow closure type so they coexist with the soil (3D, XYZ) methods of the same
# generic kernel names without ambiguity.

@kernel inbounds = true function energy_to_temperature_kernel!(out, grid, fields, closure::SnowEnergyTemperatureClosure, args...)
    i, j = @index(Global, NTuple)
    energy_to_temperature!(out, i, j, grid, fields, closure, args...)
end

@kernel inbounds = true function temperature_to_energy_kernel!(out, grid, fields, closure::SnowEnergyTemperatureClosure, args...)
    i, j = @index(Global, NTuple)
    temperature_to_energy!(out, i, j, grid, fields, closure, args...)
end
