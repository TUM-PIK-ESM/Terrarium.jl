"""
    $TYPEDEF

A collection of general physical constants that do not (usually) need to be varied in parameter calibration.
"""
@kwdef struct PhysicalConstants{NF}
    "Density of water in kg/m^3"
    ρw::NF = 1000.0

    "Density of ice in kg/m^3"
    ρi::NF = 916.2

    "Sepcific latent heat of fusion of water in J/kg"
    Lsl::NF = 3.34e5
    
    "Specific latent heat of vaporization of water in J/kg"
    Lsg::NF = 2.257e6
    
    "Gravitational constant in m/s^2"
    g::NF = 9.80665

    "Reference temperature (0°C in Kelvin)"
    T0::NF = 273.15
end

PhysicalConstants(::Type{NF}; kwargs...) where {NF} = PhysicalConstants{NF}(; kwargs...)
