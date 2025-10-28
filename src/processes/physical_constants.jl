"""
    $TYPEDEF

A collection of general physical constants that do not (usually) need to be varied in parameter calibration.
"""
@kwdef struct PhysicalConstants{NF}
    "Density of water in kg/m^3"
    ρw::NF = 1000.0

    "Density of ice in kg/m^3"
    ρi::NF = 916.2

    "Density of air at standard pressure and 0°C in kg/m^3"
    ρₐ::NF = 1.293u"kg/m^3"

    "Sepcific latent heat of fusion of water in J/kg"
    Lsl::NF = 3.34e5
    
    "Specific latent heat of vaporization of water in J/kg"
    Lsg::NF = 2.257e6
    
    "Gravitational constant in m/s^2"
    g::NF = 9.80665

    "Reference temperature (0°C in Kelvin)"
    T0::NF = 273.15

    "Stefan-Boltzmann constant in J/(s*m^2*K^4)"
    σ::NF = 5.6704e-8
    
    "von Kármán constant"
    κ::NF = 0.4

    "Psychrometric constant"
    γ::NF = 0.622

    "Specific gas constant of air in J/(kg*K)"
    Rₐ::NF = 287.058

    "Volumetric heat capacity of dry air at standard pressure and 0°C in J/(m^3*K)"
    cₐ::NF = 1005.7 * ρₐ
end

PhysicalConstants(::Type{NF}; kwargs...) where {NF} = PhysicalConstants{NF}(; kwargs...)
