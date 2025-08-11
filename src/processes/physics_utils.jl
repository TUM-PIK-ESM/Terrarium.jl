"""
    $SIGNATURES

Compute partial pressure of oxygen from surface pressure in Pa.
"""
@inline function partial_pressure_O2(pres::NF) where NF
    # TODO Shouldn't this be in physical constants?
    pres_O2 = NF(0.209) * pres
    return pres_O2
end

"""
    $SIGNATURES

Compute partial pressure of CO2 from surface pressure and CO2 concentration in Pa.
"""
@inline function partial_pressure_CO2(pres::NF, conc_co2::NF) where NF
    pres_co2 = conc_co2 * NF(1e-6) * pres
    return pres_co2
end

"""
    $SIGNATURES

Computes the vapor pressure deficit from air temperature, specific humidity, and surface pressure.
"""
@inline function compute_vpd(T_air::NF, q_air::NF, pres::NF) where NF
    # Compute Saturation vapor pressure over water [Pa]
    e_sat_w = NF(6.1094e2) * exp(NF(17.625) * T_air / (NF(243.04) + T_air))

    # Convert air specific humidity to vapor pressure [Pa]
    q_to_e = q_air * pres / (NF(0.622) + NF(0.378) * q_air)

    # Compute vapor pressure deficit [Pa]
    # TODO is the max operation needed?
    vpd = max(e_sat_w - q_to_e, NF(0.1))

    return vpd
end
