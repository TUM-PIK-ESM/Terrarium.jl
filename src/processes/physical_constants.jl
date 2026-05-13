import Thermodynamics
import Thermodynamics.Parameters:
    AbstractThermodynamicsParameters
import Thermodynamics: air_density, cp_m

"""
    $TYPEDEF

A collection of general physical constants that do not (usually) need to be varied in parameter calibration.

Properties:
$FIELDS
"""
@kwdef struct PhysicalConstants{NF} <: AbstractThermodynamicsParameters{NF}
    "Density of water in kg/m^3"
    ρw::NF = 1000.0

    "Density of ice in kg/m^3"
    ρi::NF = 916.7

    "Isobaric specific heat capacity of dry air at standard pressure and 0°C in J/(m^3*K)"
    cp_d::NF = 1004.5

    "Isobaric specific heat capacity of ice at standard pressure and 0°C in J/(m^3*K)"
    cp_i::NF = 2070.0

    "Isobaric specific heat capacity of liquid water at standard pressure and 0°C in J/(m^3*K)"
    cp_l::NF = 4186.0

    "Isobaric specific heat capacity of water vapor at standard pressure and 0°C in J/(m^3*K)"
    cp_v::NF = 1846.0

    "Specific latent heat of fusion of water in J/kg at 0°C"
    Lsl::NF = 3.34e5

    "Specific latent heat of vaporization of water in J/kg at 0°C"
    Llg::NF = 2.257e6

    "Specific latent heat of sublimation of water in J/kg at 0°C"
    Lsg::NF = 2.834e6

    "Gravitational constant in m/s^2"
    g::NF = 9.80665

    "Reference temperature (0°C in Kelvin)"
    T_ref::NF = 273.16

    "Freezing temperature of water in Kelvin"
    T_freeze::NF = 273.16

    "Triple point temperature of water in Kelvin"
    T_triple::NF = 273.16

    "Triple point pressure of water in Pa"
    press_triple::NF = 611.657

    "Stefan-Boltzmann constant in J/(s*m^2*K^4)"
    σ::NF = 5.6704e-8

    "von Kármán constant"
    κ::NF = 0.4

    "Specific gas constant of dry air in J/(kg*K)"
    R_d::NF = 287.058

    "Specific gas constant of water vapor in J/(kg*K)"
    R_v::NF = 461.5

    "Atomic mass of carbon [gC/mol]"
    C_mass::NF = 12.0
end

PhysicalConstants(::Type{NF}; kwargs...) where {NF} = PhysicalConstants{NF}(; kwargs...)

@inline Thermodynamics.Parameters.R_d(c::PhysicalConstants) = c.R_d
@inline Thermodynamics.Parameters.R_v(c::PhysicalConstants) = c.R_v
@inline Thermodynamics.Parameters.cp_d(c::PhysicalConstants) = c.cp_d
@inline Thermodynamics.Parameters.cp_i(c::PhysicalConstants) = c.cp_i
@inline Thermodynamics.Parameters.cp_l(c::PhysicalConstants) = c.cp_l
@inline Thermodynamics.Parameters.cp_v(c::PhysicalConstants) = c.cp_v
@inline Thermodynamics.Parameters.LH_v0(c::PhysicalConstants) = c.Llg
@inline Thermodynamics.Parameters.LH_s0(c::PhysicalConstants) = c.Lsg
@inline Thermodynamics.Parameters.T_0(c::PhysicalConstants) = c.T_ref
@inline Thermodynamics.Parameters.T_freeze(c::PhysicalConstants) = c.T_freeze
@inline Thermodynamics.Parameters.T_triple(c::PhysicalConstants) = c.T_triple
@inline Thermodynamics.Parameters.press_triple(c::PhysicalConstants) = c.press_triple

# Derived parameters
@inline Thermodynamics.Parameters.Rv_over_Rd(c::PhysicalConstants) = Thermodynamics.Parameters.R_v(c) / Thermodynamics.Parameters.R_d(c)
@inline ε(c::PhysicalConstants) = R_d(c) / R_v(c)

"""
    specific_heat_capacity_moist_air(c::PhysicalConstants, q)

Compute the isobaric specific heat capacity [J/(kg*K)] of moist air as a function 
of the total specific humidity `q` [kg/kg]. Wrapper around 
[`cp_m`](@extref Thermodynamics.cp_m). 
"""
@inline specific_heat_capacity_moist_air(c::PhysicalConstants, q) = cp_m(c, q)
"""
    celsius_to_kelvin(c::PhysicalConstants, T)

Convert the given temperature in °C to Kelvin based on the constant `Tref`.
"""
@inline celsius_to_kelvin(c::PhysicalConstants, T) = T + c.T_ref

"""
    stefan_boltzmann(c::PhysicalConstants, T, ϵ)

Stefan-Boltzmann law ``M = \\epsilon \\sigma T^4`` where T is the surface temperature in Kelvin
and ϵ is the emissivity.
"""
@inline stefan_boltzmann(c::PhysicalConstants, T, ϵ) = ϵ * c.σ * T^4

"""
    psychrometric_constant(c::PhysicalConstants, p)

Calcualte the psychrometric constant at the given atmospheric pressure `p`.
"""
@inline psychrometric_constant(c::PhysicalConstants, p) = cp_d(c) * p / (LH_v0(c) * ε(c))
