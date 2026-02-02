"""
    $TYPEDEF

Defines the constitutive relationship between the the internal energy and temperature of a
soil volume, i.e.
```math
U(T) = T\\times C(T) - L_f \\theta_{wi} (1 - F(T))
```
where T is temperature, C(T) is the temperature-dependent heat capacity, L_f is the
volumetric latent heat of fusion, and F(T) is the constitutive relation between temperature
and the unfrozen fraction of pore water. Note that, under this formulation, zero energy corresponds to
0°C with no ice, i.e. all pore water fully thawed.

The closure relation is defined as being a mapping from the conserved quantity (energy) to the continuous
quantity (temperature), i.e. the inverse of U(T).
"""
struct SoilEnergyTemperatureClosure <: AbstractSoilEnergyClosure end

"""
Defines `temperature` as the closure variable for `SoilEnergyTemperatureClosure`.
"""
variables(::SoilEnergyTemperatureClosure) = (
    auxiliary(:temperature, XYZ(), units=u"°C", desc="Temperature of the soil volume in °C"),
    auxiliary(:liquid_water_fraction, XYZ(), domain=UnitInterval(), desc="Fraction of unfrozen water in the pore space"),
)

@propagate_inbounds function temperature_to_energy!(
    out, i, j, k, grid, fields,
    ::SoilEnergyTemperatureClosure,
    ::FreeWater,
    energy::SoilEnergyBalance{NF, OP, SoilEnergyTemperatureClosure},
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants
) where {NF, OP}
    T = fields.temperature[i, j, k] # assumed given
    L = constants.ρw*constants.Lsl
    por = porosity(i, j, k, grid, fields, strat, bgc)
    sat = saturation_water_ice(i, j, k, grid, fields, hydrology)
    # calculate unfrozen water content from temperature
    # N.B. For the free water freeze curve, the mapping from temperature to unfrozen water content
    # within the phase change region is indeterminate since it is assumed that T = 0. As such, we
    # have to assume here that the liquid water fraction is zero if T < 0 and one otherwise. This method
    # should therefore only be used for initialization and should **not** be involved in the calculation
    # of tendencies.
    liq = out.liquid_water_fraction[i, j, k] = ifelse(
        T >= zero(T),
        one(sat),
        zero(sat),
    )
    # add liquid water fraction to fields
    fields = merge(fields, (; liquid_water_fraction = out.liquid_water_fraction))
    solid = soil_matrix(i, j, k, grid, fields, strat, bgc)
    soil = SoilVolume(por, sat, liq, solid)
    C = heat_capacity(energy.thermal_properties, soil)
    # compute energy from temperature, heat capacity, and ice fraction
    U = out.internal_energy[i, j, k] = T * C - L * sat * por * (1 - liq)
    return U
end

@propagate_inbounds function energy_to_temperature!(
    out, i, j, k, grid, fields,
    ::SoilEnergyTemperatureClosure,
    fc::FreeWater,
    energy::SoilEnergyBalance,
    hydrology::AbstractSoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants
)

    U = fields.internal_energy[i, j, k] # assumed given
    L = constants.ρw*constants.Lsl
    por = porosity(i, j, k, grid, fields, strat, bgc)
    sat = saturation_water_ice(i, j, k, grid, fields, hydrology)
    Lθ = L*sat*por
    # calculate unfrozen water content
    liq = out.liquid_water_fraction[i, j, k] = liquid_water_fraction(fc, U, Lθ, sat)
    # add liquid water fraction to fields
    fields = merge(fields, (; liquid_water_fraction = out.liquid_water_fraction))
    # calculate soil volumetric fractions
    solid = soil_matrix(i, j, k, grid, fields, strat, bgc)
    soil = SoilVolume(por, sat, liq, solid)
    C = heat_capacity(energy.thermal_properties, soil)
    # calculate temperature from internal energy and liquid water fraction
    T = out.temperature[i, j, k] = energy_to_temperature(fc, U, Lθ, C)
    return T
end

"""
Calculate the unfrozen water content from the given internal energy, latent heat content, and saturation.
"""
@inline function liquid_water_fraction(::FreeWater, U::NF, Lθ::NF, sat::NF) where {NF}
    return if U >= zero(U)
        # Case 1: U ≥ Lθ -> thawed
        one(sat)
    else
        # Case 2a: -Lθ ≤ U ≤ 0 -> phase change
        # Case 2b: U < -Lθ -> frozen (zero)
        (U >= -Lθ)*(one(sat) - safediv(U, -Lθ))
    end
end

"""
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

# Kernels

@kernel inbounds=true function temperature_to_energy_kernel!(out, grid, fields, args...)
    i, j, k = @index(Global, NTuple)
    temperature_to_energy!(out, i, j, k, grid, fields, args...)
end

@kernel inbounds=true function energy_to_temperature_kernel!(out, grid, fields, args...)
    i, j, k = @index(Global, NTuple)
    energy_to_temperature!(out, i, j, k, grid, fields, args...)
end
