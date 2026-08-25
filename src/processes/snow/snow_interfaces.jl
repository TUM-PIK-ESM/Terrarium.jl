"""
    $TYPEDSIGNATURES

Launch [`compute_snow_interface_fluxes!`](@ref) to diagnose the snow↔surface/soil coupling fluxes from the
snow state and the surface energy balance outputs: the blended soil-top heat flux (`soil_heat_flux`) and
the snow surface sublimation rate (`sublimation`, the snow-fraction bulk-aerodynamic vapor flux at the
converged skin temperature). Must run after the surface energy balance (which sets `ground_heat_flux` and
the skin temperature). No-op when there is no snowpack (`snow === nothing`).
"""
compute_snow_interface_fluxes!(state, grid, ::Nothing, args...) = nothing
function compute_snow_interface_fluxes!(
        state, grid,
        snow::SingleLayerSnow,
        seb::AbstractSurfaceEnergyBalance,
        soil::AbstractSoil,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere
    )
    # atmosphere fields (air temperature/pressure/humidity, windspeed) are needed for the sublimation flux
    out = (
        sublimation = state.sublimation,
        basal_heat_flux = state.basal_heat_flux,
    )
    fields = get_fields(state, snow, seb, soil, atmos)
    launch!(grid, XY, compute_snow_interface_fluxes_kernel!, out, fields, snow, seb, soil, constants, atmos)
    return nothing
end

"""
    $TYPEDSIGNATURES

Blended soil-top heat flux [W/m²] at grid cell `i, j`: the snow-cover-fraction-weighted combination of
the snow→soil basal conductive flux `Q_base` and the bare-ground surface energy balance closure flux `G`
(`ground_heat_flux`), `f_snow·Q_base + (1 − f_snow)·G`. The bulk snow thermal conductivity entering
`Q_base` is recovered lazily from the density scheme rather than stored.
"""
@propagate_inbounds function compute_snow_soil_heat_flux(
        i, j, grid, fields,
        snow::SingleLayerSnow,
        constants::PhysicalConstants,
        soil::AbstractSoil
    )
    G = fields.ground_heat_flux[i, j]
    f = fields.snow_cover_fraction[i, j]
    d_snow = fields.snow_depth[i, j]
    T_snow = fields.snow_temperature[i, j]
    T_soil = fields.ground_temperature[i, j]
    ρ_snow = compute_snow_density(i, j, grid, fields, snow.density)
    κ_snow = compute_thermal_conductivity(snow, constants.material, ρ_snow)
    Q_base = compute_snow_basal_heat_flux(κ_snow, T_soil, T_snow, d_snow, snow.min_conduction_thickness)
    return f * Q_base + (one(f) - f) * G
end

"""
    $TYPEDSIGNATURES

Diagnose the snow↔surface/soil coupling fluxes at grid cell `i, j`, run after the surface energy balance:
the blended soil-top heat flux (see [`compute_snow_soil_heat_flux`](@ref)) and the snow surface
sublimation rate (see [`compute_snow_sublimation_flux`](@ref), the same snow-fraction vapor flux the
surface energy balance uses for the latent-flux partition).
"""
@propagate_inbounds function compute_snow_interface_fluxes!(
        out, i, j, grid, fields,
        snow::SingleLayerSnow,
        seb::AbstractSurfaceEnergyBalance,
        soil::AbstractSoil,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere
    )
    f = snow_cover_fraction(i, j, grid, fields, snow)
    E_subl = compute_snow_sublimation_flux(i, j, grid, fields, snow, atmos, constants, seb.skin_temperature)
    out.sublimation[i, j, 1] = f * E_subl
    out.basal_heat_flux[i, j, 1] = compute_snow_soil_heat_flux(i, j, grid, fields, snow, constants, soil)
    return nothing
end

@kernel inbounds = true function compute_snow_interface_fluxes_kernel!(out, grid, fields, snow::SingleLayerSnow, args...)
    i, j = @index(Global, NTuple)
    compute_snow_interface_fluxes!(out, i, j, grid, fields, snow, args...)
end
