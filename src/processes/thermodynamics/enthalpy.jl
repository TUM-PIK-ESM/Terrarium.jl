# Material-agnostic enthalpy-temperature relations
#
# These scalar maps relate volumetric internal energy `U` (J/m³) to the unfrozen (liquid) water
# fraction and to temperature, given the volumetric latent-heat content `Lθ` (J/m³) and the
# volumetric heat capacity `C` (J/m³/K).

"""
($TYPEDSIGNATURES)
Calculate the unfrozen water content from the given internal energy, latent heat content, and saturation.
"""
@inline function liquid_water_fraction(::FreeWater, U::NF, Lθ::NF, sat::NF) where {NF}
    # Case 1: U ≥ 0 -> thawed (liq = 1)
    # Case 2a: -Lθ ≤ U < 0 -> phase change (liq = 1 - U/(-Lθ))
    # Case 2b: U < -Lθ -> frozen (liq = 0), enforced by the (U ≥ -Lθ) factor.
    return ifelse(
        U >= zero(U),
        one(sat),
        (U >= -Lθ) * (one(sat) - safediv(U, -Lθ)),
    )
end

"""
($TYPEDSIGNATURES)
Calculate the inverse enthalpy function given the internal energy, latent heat content, and heat
capacity under the free water freezing characteristic.
"""
@inline function energy_to_temperature(::FreeWater, U::NF, Lθ::NF, C::NF) where {NF}
    # Case 1:  U < -Lθ      → frozen      (T = (U + Lθ)/C)
    # Case 2a: U ≥ 0        → thawed      (T = U/C)
    # Case 2b: -Lθ ≤ U < 0  → phase change (T = 0)
    return ifelse(
        U < -Lθ,
        (U + Lθ) / C,
        ifelse(U >= zero(U), U / C, zero(NF)),
    )
end
