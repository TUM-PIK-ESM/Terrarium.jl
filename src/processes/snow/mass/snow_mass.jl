# Snow mass balance, meltwater outflow, and the depth-integrated energy/mass tendencies.

"""
    $TYPEDSIGNATURES

Snow-surface sublimation rate [m/s SWE] at grid cell `i, j`, area-weighted by the snow-covered fraction
`f_snow`. The snow surface is treated as saturated: a bulk-aerodynamic vapor flux `Δq/rₐ` evaluated at
the skin temperature, with `Δq` taken over ice for a sub-freezing surface (the saturation humidity
already dispatches over ice for `T ≤ 0` — see [`saturation_specific_humidity_vapor`](@ref)). The
water-vapor mass flux `ρₐ·Δq/rₐ` is converted to a snow-water-equivalent rate via `ρ_w`. Zero without
snow (`snow === nothing`).
"""
@propagate_inbounds compute_snow_sublimation_flux(i, j, grid, fields, ::Nothing, atmos, constants, skinT) = zero(eltype(grid))
@propagate_inbounds function compute_snow_sublimation_flux(
        i, j, grid, fields,
        snow::AbstractSnow,
        atmos::AbstractAtmosphere,
        constants::PhysicalConstants,
        skinT::AbstractSkinTemperature
    )
    f = snow_cover_fraction(i, j, grid, fields, snow)
    Tₛ = skin_temperature(i, j, grid, fields, skinT)
    rₐ = aerodynamic_resistance(i, j, grid, fields, atmos)
    Δq = compute_specific_humidity_difference(i, j, grid, fields, atmos, constants, Tₛ) # over ice for Tₛ ≤ 0
    Tₐ = air_temperature(i, j, grid, fields, atmos)
    pres = air_pressure(i, j, grid, fields, atmos)
    q_air = specific_humidity(i, j, grid, fields, atmos)
    ρₐ = Thermodynamics.air_density(constants.thermodynamics, celsius_to_kelvin(constants.thermodynamics, Tₐ), pres, q_air)
    ρ_w = constants.material.density_water
    # saturated vapor mass flux converted to a snow-water-equivalent rate, area-weighted by the snow-covered
    # fraction `f` to give the grid-cell-mean sublimation (W_snow and Ū_snow are grid-cell means)
    return f * ρₐ * (Δq / rₐ) / ρ_w
end

# Kernel functions

"""
    $TYPEDSIGNATURES

Snow water equivalent (SWE) tendency (m/s) at grid cell `i, j`:
```
dW_snow/dt = S + R_snow − M − E_subl
```
where `S` is snowfall, `R_snow = f_snow · rainfall` the rain intercepted by the snow-covered fraction,
`M` the Darcy meltwater outflow (see [`snow_meltwater_flux`](@ref)), and `E_subl` the sublimation rate.
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
    E_subl = fields.sublimation[i, j]
    dWdt = S + R_snow - M - E_subl
    return dWdt
end
