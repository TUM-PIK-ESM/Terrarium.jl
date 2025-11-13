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

"""
    relative_to_specific_humidity(r_h, pr, Ts, Tair)

Derives specific humidity from measured relative humidity, air pressure, and soil/air temperatures.
"""
@inline relative_to_specific_humidity(r_h, pr, Ts, Tair) = 0.622 * (r_h / 100) * saturation_vapor_pressure(Tair, Ts) / pr

# saturation vapor pressure
"""
    saturation_vapor_pressure(T, a₁, a₂, a₃)

August-Roche-Magnus equation for saturation vapor pressure at temperature `T` with empirical
coefficients a₁, a₂, and a₃.
"""
@inline saturation_vapor_pressure(T, a₁, a₂, a₃) = a₁ * exp(a₂ * T / (T + a₃))

"""
    saturation_vapor_pressure(T, Ts=T)

Saturation vapor pressure at the given temperature `T` over a surface at temperature `Ts`,
accounting for both frozen (`Ts < 0°C`) and unfrozen conditions.

Coefficients taken from Alduchov and Eskridge (1997).
"""
@inline function saturation_vapor_pressure(T::NF, Ts::NF=T) where {NF}
    if Ts <= zero(Ts)
        saturation_vapor_pressure(T, NF(611.0), NF(22.46), NF(272.62))
    else
        saturation_vapor_pressure(T, NF(611.0), NF(17.62), NF(243.12))
    end
end
