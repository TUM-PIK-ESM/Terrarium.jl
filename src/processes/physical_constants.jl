import Thermodynamics
import Thermodynamics.Parameters:
    AbstractThermodynamicsParameters
import Thermodynamics: air_density, cp_m

@kwdef struct ThermodynamicConstants{NF} <: AbstractThermodynamicsParameters{NF}
    "Isobaric specific heat capacity of dry air at standard pressure and 0°C in J/(m^3*K)"
    specific_heat_capacity_dry_air::NF = 1004.5
    "Isobaric specific heat capacity of ice at standard pressure and 0°C in J/(m^3*K)"
    specific_heat_capacity_ice::NF = 2070.0
    "Isobaric specific heat capacity of liquid water at standard pressure and 0°C in J/(m^3*K)"
    specific_heat_capacity_liquid_water::NF = 4186.0
    "Isobaric specific heat capacity of water vapor at standard pressure and 0°C in J/(m^3*K)"
    specific_heat_capacity_water_vapor::NF = 1846.0
    "Specific latent heat of fusion of water in J/kg at 0°C"
    latent_heat_fusion_at_reference::NF = 3.34e5
    "Specific latent heat of vaporization of water in J/kg at 0°C"
    latent_heat_vaporization_at_reference::NF = 2.257e6
    "Specific latent heat of sublimation of water in J/kg at 0°C"
    latent_heat_sublimation_at_reference::NF = 2.834e6
    "Reference temperature (0°C in Kelvin)"
    temperature_reference::NF = 273.16
    "Freezing temperature of water in Kelvin"
    temperature_water_freeze::NF = 273.16
    "Triple point temperature of water in Kelvin"
    temperature_water_triple_point::NF = 273.16
    "Triple point pressure of water in Pa"
    pressure_water_triple_point::NF = 611.657
    "Specific gas constant of dry air in J/(kg*K)"
    gas_constant_dry_air::NF = 287.058
    "Specific gas constant of water vapor in J/(kg*K)"
    gas_constant_water_vapor::NF = 461.5
end

ThermodynamicConstants(::Type{NF}; kwargs...) where {NF} = ThermodynamicConstants{NF}(; kwargs...)


@kwdef struct MaterialConstants{NF}
    "Density of water in kg/m^3"
    water_density::NF = 1000.0
    "Density of ice in kg/m^3"
    ice_density::NF = 916.7
    "Atomic mass of carbon [gC/mol]"
    atomic_weight_carbon::NF = 12.0
end

MaterialConstants(::Type{NF}; kwargs...) where {NF} = MaterialConstants{NF}(; kwargs...)

@kwdef struct UniversalConstants{NF}
    "Gravitational constant in m/s^2"
    gravitational_acceleration::NF = 9.80665
    "Stefan-Boltzmann constant in J/(s*m^2*K^4)"
    stefan_boltzmann_constant::NF = 5.6704e-8
    "von Kármán constant"
    von_karman_constant::NF = 0.4
end

UniversalConstants(::Type{NF}; kwargs...) where {NF} = UniversalConstants{NF}(; kwargs...)

"""
    $TYPEDEF

A collection of general physical constants that do not (usually) need to be varied in parameter calibration.

Properties:
$FIELDS
"""
struct PhysicalConstants{NF}
    thermodynamic_constants::ThermodynamicConstants{NF}
    material_constants::MaterialConstants{NF}
    universal_constants::UniversalConstants{NF}
end

PhysicalConstants() = PhysicalConstants(Float64)

function PhysicalConstants(
        ::Type{NF};
        thermodynamic_constants::ThermodynamicConstants{NF} = ThermodynamicConstants(NF),
        material_constants::MaterialConstants{NF} = MaterialConstants(NF),
        universal_constants::UniversalConstants{NF} = UniversalConstants(NF),
    ) where {NF}
    return PhysicalConstants{NF}(thermodynamic_constants, material_constants, universal_constants)
end

@inline Thermodynamics.Parameters.R_d(c::ThermodynamicConstants) = c.gas_constant_dry_air
@inline Thermodynamics.Parameters.R_v(c::ThermodynamicConstants) = c.gas_constant_water_vapor
@inline Thermodynamics.Parameters.cp_d(c::ThermodynamicConstants) = c.specific_heat_capacity_dry_air
@inline Thermodynamics.Parameters.cp_i(c::ThermodynamicConstants) = c.specific_heat_capacity_ice
@inline Thermodynamics.Parameters.cp_l(c::ThermodynamicConstants) = c.specific_heat_capacity_liquid_water
@inline Thermodynamics.Parameters.cp_v(c::ThermodynamicConstants) = c.specific_heat_capacity_water_vapor
@inline Thermodynamics.Parameters.LH_v0(c::ThermodynamicConstants) = c.latent_heat_vaporization_at_reference
@inline Thermodynamics.Parameters.LH_s0(c::ThermodynamicConstants) = c.latent_heat_sublimation_at_reference
@inline Thermodynamics.Parameters.T_0(c::ThermodynamicConstants) = c.temperature_reference
@inline Thermodynamics.Parameters.T_freeze(c::ThermodynamicConstants) = c.temperature_water_freeze
@inline Thermodynamics.Parameters.T_triple(c::ThermodynamicConstants) = c.temperature_water_triple_point
@inline Thermodynamics.Parameters.press_triple(c::ThermodynamicConstants) = c.pressure_water_triple_point

# Derived parameters
@inline Thermodynamics.Parameters.Rv_over_Rd(c::ThermodynamicConstants) = Thermodynamics.Parameters.R_v(c) / Thermodynamics.Parameters.R_d(c)
@inline ε(c::ThermodynamicConstants) = R_d(c) / R_v(c)

"""
    specific_heat_capacity_moist_air(c::ThermodynamicConstants, q)

Compute the isobaric specific heat capacity [J/(kg*K)] of moist air as a function 
of the total specific humidity `q` [kg/kg]. Wrapper around 
[`cp_m`](@extref Thermodynamics.cp_m). 
"""
@inline specific_heat_capacity_moist_air(c::ThermodynamicConstants, q) = cp_m(c, q)
"""
    celsius_to_kelvin(c::ThermodynamicConstants, T)

Convert the given temperature in °C to Kelvin based on the constant `temperature_reference`.
"""
@inline celsius_to_kelvin(c::ThermodynamicConstants, T) = T + c.temperature_reference

"""
    stefan_boltzmann(c::UniversalConstants, T, ϵ)

Stefan-Boltzmann law ``M = \\epsilon \\sigma T^4`` where T is the surface temperature in Kelvin
and ϵ is the emissivity and σ is the Stefan-Boltzmann constant.
"""
@inline stefan_boltzmann(c::UniversalConstants, T, ϵ) = ϵ * c.stefan_boltzmann_constant * T^4

"""
    psychrometric_constant(c::ThermodynamicConstants, p)

Calcualte the psychrometric constant at the given atmospheric pressure `p`.
"""
@inline psychrometric_constant(c::ThermodynamicConstants, p) = c.specific_heat_capacity_dry_air * p / (c.latent_heat_vaporization_at_reference * ε(c))
