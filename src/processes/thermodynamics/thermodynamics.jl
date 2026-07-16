# Thermodynamics interface

@inline Thermodynamics.Parameters.R_d(c::ThermodynamicConstants) = c.gas_constant_dry_air
@inline Thermodynamics.Parameters.R_v(c::ThermodynamicConstants) = c.gas_constant_water_vapor
@inline Thermodynamics.Parameters.cp_d(c::ThermodynamicConstants) = c.specific_heat_capacity_dry_air
@inline Thermodynamics.Parameters.cp_i(c::ThermodynamicConstants) = c.specific_heat_capacity_ice
@inline Thermodynamics.Parameters.cp_l(c::ThermodynamicConstants) = c.specific_heat_capacity_liquid_water
@inline Thermodynamics.Parameters.cp_v(c::ThermodynamicConstants) = c.specific_heat_capacity_water_vapor
@inline Thermodynamics.Parameters.LH_v0(c::ThermodynamicConstants) = c.latent_heat_vaporization
@inline Thermodynamics.Parameters.LH_s0(c::ThermodynamicConstants) = c.latent_heat_sublimation
@inline Thermodynamics.Parameters.T_0(c::ThermodynamicConstants) = c.temperature_reference
@inline Thermodynamics.Parameters.T_freeze(c::ThermodynamicConstants) = c.temperature_water_freeze
@inline Thermodynamics.Parameters.T_triple(c::ThermodynamicConstants) = c.temperature_water_triple_point
@inline Thermodynamics.Parameters.press_triple(c::ThermodynamicConstants) = c.pressure_water_triple_point

# Derived parameters
@inline Thermodynamics.Parameters.Rv_over_Rd(c::ThermodynamicConstants) = Thermodynamics.Parameters.R_v(c) / Thermodynamics.Parameters.R_d(c)
@inline ratio_gas_constant_dry_air_to_water_vapor(c::ThermodynamicConstants) = R_d(c) / R_v(c)

"""
    specific_heat_capacity_moist_air(c::ThermodynamicConstants, q)

Compute the isobaric specific heat capacity [J/(kg*K)] of moist air as a function 
of the total specific humidity `q` [kg/kg]. Wrapper around 
[`cp_m`](@extref Thermodynamics.cp_m). 
"""
@inline specific_heat_capacity_moist_air(c::ThermodynamicConstants, q) = Thermodynamics.cp_m(c, q)

"""
    celsius_to_kelvin(c::ThermodynamicConstants, T)

Convert the given temperature in °C to Kelvin based on the constant `temperature_water_freeze`.
"""
@inline celsius_to_kelvin(c::ThermodynamicConstants, T) = T + c.temperature_water_freeze

"""
    psychrometric_constant(c::ThermodynamicConstants, p)

Calcualte the psychrometric constant at the given atmospheric pressure `p`.
"""
@inline psychrometric_constant(c::ThermodynamicConstants, p) = c.specific_heat_capacity_dry_air * p / (c.latent_heat_vaporization * ratio_gas_constant_dry_air_to_water_vapor(c))

"""
    relative_to_specific_humidity(r_h, pr, T, c::ThermodynamicConstants)

Derives specific humidity from measured relative humidity `r_h` [%], air pressure `pr` [Pa],
air temperature `T` [°C], and physical constants `c`. Assumes saturation over ice for
`T <= 0°C` and over liquid water otherwise. Wrapper around
[`q_vap_from_RH`](@extref Thermodynamics.q_vap_from_RH).
"""
@inline function relative_to_specific_humidity(r_h, pr, T, c::ThermodynamicConstants)
    T_K = celsius_to_kelvin(c, T)
    phase = T <= zero(T) ? Thermodynamics.Ice() : Thermodynamics.Liquid()
    return Thermodynamics.q_vap_from_RH(c, pr, T_K, r_h / 100, phase)
end

"""
    vapor_pressure_to_specific_humidity(c::ThermodynamicConstants, e, pr)

Derives specific humidity from measured vapor pressure `e` [Pa] and air pressure `pr` [Pa]. 
"""
@inline function vapor_pressure_to_specific_humidity(c::ThermodynamicConstants, e, pr)
    return ratio_gas_constant_dry_air_to_water_vapor(c) * e / (pr - e * (1 - ratio_gas_constant_dry_air_to_water_vapor(c)))
end

"""
    saturation_vapor_pressure(T)

Saturation vapor pressure of an air parcel at the given temperature `T` in °C. By default, the saturation vapor
pressure is computed over ice for `T <= 0°C` and over water for `T > 0°C`. Wrapper around
[`saturation_vapor_pressure`](@extref Thermodynamics.saturation_vapor_pressure).
"""
@inline function saturation_vapor_pressure(c::ThermodynamicConstants, T::NF) where {NF}
    T_K = celsius_to_kelvin(c, T)
    return ifelse(
        T <= zero(T),
        Thermodynamics.saturation_vapor_pressure(c, T_K, Thermodynamics.Ice()),
        Thermodynamics.saturation_vapor_pressure(c, T_K, Thermodynamics.Liquid())
    )
end

"""
    saturation_specific_humidity_vapor(c::ThermodynamicConstants, T, ρ)

Saturation specific humidity at temperature `T` [°C] and density `ρ` [kg/m³]. Dispatches
over ice for `T <= 0°C` and over liquid water otherwise. Wrapper around
[`q_vap_saturation`](@extref Thermodynamics.q_vap_saturation).
"""
@inline function saturation_specific_humidity_vapor(c::ThermodynamicConstants, T, ρ)
    T_K = celsius_to_kelvin(c, T)
    return ifelse(
        T <= zero(T),
        Thermodynamics.q_vap_saturation(c, T_K, ρ, Thermodynamics.Ice()),
        Thermodynamics.q_vap_saturation(c, T_K, ρ, Thermodynamics.Liquid())
    )
end

"""
    stefan_boltzmann(c::UniversalConstants, T, ϵ)

Stefan-Boltzmann law ``M = \\epsilon \\sigma T^4`` where T is the surface temperature in Kelvin
and ϵ is the emissivity and σ is the Stefan-Boltzmann constant.
"""
@inline stefan_boltzmann(c::UniversalConstants, T, ϵ) = ϵ * c.stefan_boltzmann_constant * T^4
