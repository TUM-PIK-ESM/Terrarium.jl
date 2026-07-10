"""
Return the number of seconds per day in the given number format.
"""
seconds_per_day(::Type{NF}) where {NF} = ustrip(u"s", NF(1)u"d")

"""
Return the number of seconds per hour in the given number format.
"""
seconds_per_hour(::Type{NF}) where {NF} = ustrip(u"s", NF(1)u"hr")

"""
    $SIGNATURES

Compute partial pressure of oxygen from surface pressure in Pa.
"""
@inline function partial_pressure_O2(pres::NF) where {NF}
    # TODO Shouldn't this be in physical constants?
    pres_O2 = NF(0.209) * pres
    return pres_O2
end

"""
    $SIGNATURES

Compute partial pressure of CO2 from surface pressure and CO2 concentration in Pa.
"""
@inline function partial_pressure_CO2(pres::NF, conc_co2::NF) where {NF}
    pres_co2 = conc_co2 * NF(1.0e-6) * pres
    return pres_co2
end
