@kwdef struct PhyscialConstants{NF}
    "Density of water in kg/m^3"
    ρw::NF = 1000.0

    "Density of ice in kg/m^3"
    ρi::NF = 916.2

    "Sepcific latent heat of fusion of water in J/kg"
    # Lsl = 3.34e5u"J/kg", # specific latent heat of fusion of water [J/kg]
    # Lsg = 2.257e6u"J/kg", # specific latent heat of vaporization of water [J/kg]
    # g = 9.80665u"m/s^2", # gravitational constant  
end
