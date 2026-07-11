# Snow geometric and thermal properties for the simple single-layer scheme.
#
# These are pure functions of the snow water equivalent `W` and the constant bulk density `ρ_s`; for a
# fixed `ρ_s` the thermal conductivity is constant. They are the "properties" layer of the scheme and
# carry no coupling to other processes.

"""
    $TYPEDSIGNATURES

Return the bulk snow density `ρ_s` [kg/m³] from the snow process's density scheme.
"""
@inline snow_density(snow::SingleLayerSnow) = snow_density(snow.density)

"""
    $TYPEDSIGNATURES

Snow layer depth `d_s = W·ρ_w/ρ_s` [m], converting the water-equivalent depth `W` [m] to the physical
snow depth using the water density `ρ_w` and the bulk snow density `ρ_s`.
"""
@inline compute_snow_depth(snow::SingleLayerSnow, W::NF, ρ_w::NF) where {NF} = max(W, zero(NF)) * ρ_w / snow_density(snow)

"""
    $TYPEDSIGNATURES

Sub-grid snow-covered area fraction `f_snow = W/(W + W_ref)` ∈ [0,1). Smooth and differentiable, with
`f_snow → 0` as `W → 0` and `f_snow → 1` as `W → ∞`.
"""
@inline function compute_snow_cover_fraction(snow::SingleLayerSnow, W::NF) where {NF}
    Wp = max(W, zero(NF))
    return Wp / (Wp + snow.cover_reference)
end

"""
    $TYPEDSIGNATURES

Bulk snow thermal conductivity via the density power law `κ_snow = a·(ρ_s/ρ_w)^b`
([yen1981review](@cite)), constant for a fixed bulk density `ρ_s`.
"""
@inline compute_snow_thermal_conductivity(snow::SingleLayerSnow, ρ_w::NF) where {NF} =
    snow.conductivity_coefficient * (snow_density(snow) / ρ_w)^snow.conductivity_exponent

"""
    $TYPEDSIGNATURES

Volumetric snow internal energy `U_v = E/d_s` [J/m³] from the depth-integrated energy `E` [J/m²] and the
snow depth `d_s` [m]. Guarded with [`safediv`](@ref); the indeterminate `W → 0` limit is masked
downstream by `f_snow → 0`.
"""
@inline volumetric_snow_energy(E::NF, d_s::NF) where {NF} = safediv(E, d_s)

# Field-level accessors

@propagate_inbounds snow_water_equivalent(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_water_equivalent[i, j]

@propagate_inbounds snow_energy(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_energy[i, j]

@propagate_inbounds snow_depth(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_depth[i, j]

@propagate_inbounds snow_cover_fraction(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_cover_fraction[i, j]

@propagate_inbounds snow_thermal_conductivity(i, j, grid, fields, ::SingleLayerSnow) = fields.snow_thermal_conductivity[i, j]

# NoSnow: inert — no snowpack present, so all field-level accessors return zero.

@propagate_inbounds snow_water_equivalent(i, j, grid, fields, ::NoSnow) = zero(eltype(grid))

@propagate_inbounds snow_energy(i, j, grid, fields, ::NoSnow) = zero(eltype(grid))

@propagate_inbounds snow_depth(i, j, grid, fields, ::NoSnow) = zero(eltype(grid))

@propagate_inbounds snow_cover_fraction(i, j, grid, fields, ::NoSnow) = zero(eltype(grid))

@propagate_inbounds snow_thermal_conductivity(i, j, grid, fields, ::NoSnow) = zero(eltype(grid))
