# Material-agnostic enthalpy-temperature relations
#
# These scalar maps relate volumetric internal energy `U` (J/m³) to the unfrozen (liquid) water
# fraction and to temperature, given the volumetric latent-heat content `Lθ` (J/m³) and the
# volumetric heat capacity `C` (J/m³/K).

"""
    $TYPEDSIGNATURES

Calculate the unfrozen water content from the given internal energy, latent heat content, and saturation.
"""
@inline function liquid_water_fraction(::FreeWater, U::NF, Lθ::NF, sat::NF) where {NF}
    return if U >= zero(U)
        # Case 1: U ≥ Lθ -> thawed
        one(sat)
    else
        # Case 2a: -Lθ ≤ U ≤ 0 -> phase change
        # Case 2b: U < -Lθ -> frozen (zero)
        (U >= -Lθ) * (one(sat) - safediv(U, -Lθ))
    end
end

"""
    $TYPEDSIGNATURES
    
Calculate the inverse enthalpy function given the internal energy, latent heat content, and heat
capacity under the free water freezing characteristic.
"""
@inline function energy_to_temperature(::FreeWater, U::NF, Lθ::NF, C::NF) where {NF}
    return if U < -Lθ
        # Case 1: U < -Lθ → frozen
        (U + Lθ) / C
    elseif U >= zero(U)
        # Case 2a: U ≥ 0 → thawed
        U / C
    else
        # Case 2b: -Lθ ≤ U < 0 → phase change
        zero(NF)
    end
    # One-liner version:
    # return (U < -Lθ)*(U + Lθ) / C
end
