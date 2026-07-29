# Snow mass balance, meltwater outflow, and the depth-integrated energy/mass tendencies.

"""
    $TYPEDSIGNATURES

Darcy-type meltwater outflow `M_r` (m/s, SWE). Liquid water in excess of the capillary retention `L_c`
drains from the snowpack with a cubic conductivity [tarbotonSpatiallyDistributedEnergy1994, Eq. (23-24)](@cite) (in excess-saturation
form): `M_r = K_sat · S*³` with `S* = max(θ_liq − L_c, 0) / (1 − L_c)`, where `θ_liq` is the liquid
fraction of the water substance. Outflow vanishes smoothly as `θ_liq → L_c` and saturates at `K_sat` as
`θ_liq → 1`.
# References

* [tarbotonSpatiallyDistributedEnergy1994](@cite) Tarboton et al., Report (1994)
"""
@inline function compute_meltwater_outflow(hydraulics::ConstantSnowHydraulics, θ_liq::NF) where {NF}
    liq_c = hydraulics.capillary_retention
    K_sat = hydraulics.saturated_conductivity
    Sstar = max(θ_liq - liq_c, zero(NF)) / (one(NF) - liq_c)
    return K_sat * Sstar^3
end

# Kernel functions

"""
    $TYPEDSIGNATURES

Snow water equivalent (SWE) tendency [m/s] at grid cell `i, j`:
```
dW_snow/dt = S + R_snow − M − E_sub
```
where `S` is snowfall, `R_snow = f_snow · rainfall` the rain intercepted by the snow-covered fraction,
`M` the Darcy meltwater outflow (see [`snow_meltwater_flux`](@ref)), and `E_sub` the sublimation rate.
"""
@propagate_inbounds function compute_snow_water_tendency(
        i, j, grid, fields,
        snow::SingleLayerSnow,
        atmos::AbstractAtmosphere
    )
    f_snow = snow_cover_fraction(i, j, grid, fields, snow)
    S = snowfall(i, j, grid, fields, atmos)
    M = snow_meltwater_flux(i, j, grid, fields, snow)
    R = rainfall(i, j, grid, fields, atmos)
    R_snow = f_snow * R
    E_sub = fields.sublimation[i, j]
    dWdt = S + R_snow - M - E_sub
    return dWdt
end
