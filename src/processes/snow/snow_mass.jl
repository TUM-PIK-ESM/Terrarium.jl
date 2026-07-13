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

Advected heat flux [W/m²] carried into the snowpack by precipitation (UEB eqn 13), relative to ice at
0 °C: snowfall `P_s` arrives as ice at `min(T_air, 0)`, while rain-on-snow `R_on_snow` arrives as liquid
carrying its latent heat of fusion plus sensible heat at `max(T_air, 0)`.
"""
@inline function compute_rain_heat_flux(P_s::NF, R_on_snow::NF, T_air::NF, constants::PhysicalConstants) where {NF}
    ρ_w = constants.material.density_water
    L_f = constants.thermodynamics.latent_heat_fusion
    c_i = constants.thermodynamics.specific_heat_capacity_ice
    c_w = constants.thermodynamics.specific_heat_capacity_liquid_water
    Q_snow = ρ_w * P_s * c_i * min(T_air, zero(NF))
    Q_rain = ρ_w * R_on_snow * (L_f + c_w * max(T_air, zero(NF)))
    return Q_snow + Q_rain
end

# Kernel functions

"""
    $TYPEDSIGNATURES

Accumulate the snow water-equivalent and depth-integrated energy tendencies at grid cell `i, j`.

Mass balance (SWE, continuous):
```
dW/dt = P_s + R_on_snow − M_r − E_subl
```
Energy balance (depth-integrated `E`, all fluxes positive upward):
```
dE/dt = G_base − G_top + Q_precip − Q_melt
```
where `P_s` is snowfall, `R_on_snow = f_snow·rainfall` the rain intercepted by the snow-covered
fraction, `M_r` the Darcy meltwater outflow, `E_subl` the sublimation rate, `G_top`/`G_base` the
surface/basal heat fluxes, `Q_precip` the advected precipitation heat, and `Q_melt = ρ_w·L_f·M_r` the
heat leaving with meltwater at 0 °C.
"""
@propagate_inbounds function compute_snow_tendencies!(
        tendencies, i, j, grid, fields,
        snow::SingleLayerSnow,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants
    )
    θ_liq = fields.snow_liquid_fraction[i, j]
    f_snow = fields.snow_cover_fraction[i, j]
    P_s = snowfall(i, j, grid, fields, atmos)
    rain = rainfall(i, j, grid, fields, atmos)
    T_air = air_temperature(i, j, grid, fields, atmos)
    G_top = fields.surface_heat_flux[i, j]   # positive upward: energy leaving the snow top
    G_base = fields.basal_heat_flux[i, j]    # positive upward: energy entering the snow base from soil
    E_subl = fields.sublimation[i, j]
    ρ_w = constants.material.density_water
    L_f = constants.thermodynamics.latent_heat_fusion
    # rainfall intercepted by the snow-covered fraction adds to the pack; the rest bypasses to the soil
    R_on_snow = f_snow * rain
    M_r = compute_meltwater_outflow(snow.hydraulic_properties, θ_liq)
    Q_precip = compute_rain_heat_flux(P_s, R_on_snow, T_air, constants)
    Q_melt = ρ_w * L_f * M_r
    tendencies.snow_water_equivalent[i, j, 1] += P_s + R_on_snow - M_r - E_subl
    tendencies.snow_energy[i, j, 1] += G_base - G_top + Q_precip - Q_melt
    return nothing
end
