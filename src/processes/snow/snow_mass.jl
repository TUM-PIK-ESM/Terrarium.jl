# Snow mass balance, meltwater outflow, and the depth-integrated energy/mass tendencies.

"""
    $TYPEDSIGNATURES

Darcy-type meltwater outflow `M_r` [m/s] (SWE). Liquid water in excess of the capillary retention `L_c`
drains from the pack with a cubic conductivity (Male & Gray 1981; UEB eqns 23–24, in excess-saturation
form): `M_r = K_sat · S*³` with `S* = max(θ_liq − L_c, 0) / (1 − L_c)`, where `θ_liq` is the liquid
fraction of the water substance. Outflow vanishes smoothly as `θ_liq → L_c` and saturates at `K_sat` as
`θ_liq → 1`.
"""
@inline function compute_meltwater_outflow(hydraulics::ConstantSnowHydraulics, θ_liq::NF) where {NF}
    L_c = hydraulics.capillary_retention
    K_sat = hydraulics.saturated_conductivity
    Sstar = max(θ_liq - L_c, zero(NF)) / (one(NF) - L_c)
    return K_sat * Sstar^3
end

"""
    $TYPEDSIGNATURES

Advected heat flux [W/m²] carried into the snowpack by precipitation, relative to liquid water at 0 °C
(the `U = 0` reference of the `FreeWater` enthalpy closure). Fresh snow `P_s` arrives as ice, which sits
`L_f` below the liquid reference, plus sensible heat for `T_air < 0`; rain-on-snow `R_on_snow` arrives as
liquid carrying only its sensible heat for `T_air > 0`. The latent heat released when rain refreezes in a
cold pack is captured implicitly by the enthalpy closure, so it is *not* added here (adding `L_f` would
double-count relative to the liquid-water reference).
"""
@inline function compute_precip_heat_flux(P_s::NF, R_on_snow::NF, T_air::NF, constants::PhysicalConstants) where {NF}
    ρ_w = constants.material.density_water
    L_f = constants.thermodynamics.latent_heat_fusion
    c_i = constants.thermodynamics.specific_heat_capacity_ice
    c_w = constants.thermodynamics.specific_heat_capacity_liquid_water
    Q_snow = ρ_w * P_s * (c_i * min(T_air, zero(NF)) - L_f)
    Q_rain = ρ_w * R_on_snow * c_w * max(T_air, zero(NF))
    Q_prcp = Q_snow + Q_rain
    return Q_prcp
end

# Kernel functions

"""
    $TYPEDSIGNATURES

Snow water-equivalent tendency [m/s SWE] at grid cell `i, j` (continuous mass balance):
```
dW/dt = P_s + R_on_snow − M_r − E_subl
```
where `P_s` is snowfall, `R_on_snow = f_snow·rainfall` the rain intercepted by the snow-covered fraction,
`M_r` the Darcy meltwater outflow (see [`snow_meltwater_flux`](@ref)), and `E_subl` the sublimation rate.
"""
@propagate_inbounds function compute_snow_water_tendency(
        i, j, grid, fields,
        snow::SingleLayerSnow,
        atmos::AbstractAtmosphere
    )
    f_snow = snow_cover_fraction(i, j, grid, fields, snow)
    S = snowfall(i, j, grid, fields, atmos)
    M = snow_meltwater_flux(i, j, grid, fields, snow)
    R = f_snow * rainfall(i, j, grid, fields, atmos)
    E_s = fields.sublimation[i, j]
    dWdt = S + R - M - E_s
    return dWdt
end

"""
    $TYPEDSIGNATURES

Depth-integrated snow energy tendency [W/m²] at grid cell `i, j` (all fluxes positive upward):
```
dE/dt = G_base − G_top + Q_precip
```
where `G_top`/`G_base` are the surface/basal heat fluxes and `Q_precip` the advected precipitation heat
(see [`compute_precip_heat_flux`](@ref)). No explicit meltwater energy term appears: meltwater drains as
liquid water at 0 °C, which is the zero-enthalpy reference (`U = 0`) of the `FreeWater` closure, so it
carries no enthalpy out of the pack.
"""
@propagate_inbounds function compute_snow_energy_tendency(
        i, j, grid, fields,
        snow::SingleLayerSnow,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants
    )
    G_top = fields.surface_heat_flux[i, j]   # positive upward: energy leaving the snow top
    G_base = fields.basal_heat_flux[i, j]    # positive upward: energy entering the snow base from soil
    f_snow = snow_cover_fraction(i, j, grid, fields, snow)
    P_s = snowfall(i, j, grid, fields, atmos)
    R_on_snow = f_snow * rainfall(i, j, grid, fields, atmos)
    T_air = air_temperature(i, j, grid, fields, atmos)
    Q_precip = compute_precip_heat_flux(P_s, R_on_snow, T_air, constants)
    return G_base - G_top + Q_precip
end

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
