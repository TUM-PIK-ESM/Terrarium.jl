# Note: this implementaion assumes a daily timestep but this should be changed later to allow for a more flexible timestep 
"""
    $TYPEDEF
Photosynthesis implementation from PALADYN (Willeit 2016) for C3 PFTs following the general light
use efficiency model described in Haxeltine and Prentice 1996.

Authors: Maha Badri and Matteo Willeit

Properties:
$TYPEDFIELDS
"""
@kwdef struct LUEPhotosynthesis{NF} <: AbstractPhotosynthesis
    # TODO check physical meaning of this parameter + add unit
    "Value of τ at 25°C"
    τ25::NF = 2600.0 

    # TODO check physical meaning of this parameter + add unit
    "Value of Kc at 25°C"
    Kc25::NF = 30.0

    # TODO check physical meaning of this parameter + add unit
    "Value of Ko at 25°C"
    Ko25::NF = 3.0e4

    # TODO check physical meaning of this parameter + add unit
    "q10 for temperature-sensitive parameter τ"
    q10_τ::NF = 0.57

    # TODO check physical meaning of this parameter + add unit
    "q10 for temperature-sensitive parameter Kc"
    q10_Kc::NF = 2.1

    # TODO check physical meaning of this parameter + add unit
    "q10 for temperature-sensitive parameter Ko"
    q10_Ko::NF = 1.2

    "Leaf albedo in PAR range [-]"
    α_leaf::NF = 0.17

    "Conversion factor for solar radiation at 550 nm from J/m² to mol/m² [mol/J]"
    cq::NF = 4.6e-6

    "Extinction coefficient for radiation through vegetation [-]"
    k_ext::NF = 0.5

    "Fraction of PAR assimilated at ecosystem level, relative to leaf level [-]"
    α_a::NF = 0.5

    # TODO check physical meaning of this parameter 
    "Parameter, PFT specific [°C]"
    t_CO2_high::NF = 42.0 # Value for Needleleaf tree PFT 

    # TODO check physical meaning of this parameter
    "Parameter, PFT specific [°C]"
    t_CO2_low::NF = -4.0 # Value for Needleleaf tree PFT 

    # TODO check physical meaning of this parameter
    "Parameter, PFT specific [°C]"
    t_photos_high::NF = 30.0 # Value for Needleleaf tree PFT 

    # TODO check physical meaning of this parameter
    "Parameter, PFT specific [°C]"
    t_photos_low::NF = 15.0 # Value for Needleleaf tree PFT 

    "Intrinsic quantum efficiency of CO2 uptake in C3 plants [mol/mol]"
    α_C3::NF = 0.08 

    "Atomic mass of carbon [gC/mol]"
    C_mass::NF = 12.0

    # TODO check physical meaning of this parameter
    "Shape parameter [-]"
    θ_r::NF = 0.7
end

LUEPhotosynthesis(::Type{NF}; kwargs...) where {NF} = LUEPhotosynthesis{NF}(; kwargs...)

variables(::LUEPhotosynthesis{NF}) where {NF} = (
    auxiliary(:An, XY()), # Daily net assimilation (photosynthesis) [gC/m²/day]
    auxiliary(:Rd, XY()), # Daily leaf respiration [gC/m²/day]
    auxiliary(:GPP, XY()), # Gross Primary Production [kgC/m²/day]
    input(:SMLF, XY(), default=NF(1)), # soil moisture limiting factor with default value of 1
    input(:λc, XY()), # Ratio of leaf-internal to air CO2 concentration [-]
    input(:LAI, XY()), # Leaf Area Index [m²/m²]
)

"""
    $SIGNATURES
Computes kinetic parameters `τ`, `Kc`, `Ko` based on temperature.
"""

@inline function compute_kinetic_parameters(photo::LUEPhotosynthesis{NF}, T_air::NF) where NF
    # TODO check meaning of these parameters, Appendix C in PALADYN paper
    # TODO add units
    τ = photo.τ25 * photo.q10_τ^((T_air - NF(25.0)) * NF(0.1))
    Kc = photo.Kc25 * photo.q10_Kc^((T_air - NF(25.0)) * NF(0.1))
    Ko = photo.Ko25 * photo.q10_Ko^((T_air - NF(25.0)) * NF(0.1))
    return τ, Kc, Ko
end


"""
    $SIGNATURES
Computes the CO2 compensation point `Γ_star`,
Eq. C6, PALADYN (Willeit 2016).
"""
@inline function compute_Γ_star(photo::LUEPhotosynthesis{NF}, τ::NF, pres_O2::NF) where NF
    # TODO add unit
    Γ_star = pres_O2 / (NF(2.0) * τ)
    return Γ_star
end

"""
    $SIGNATURES
Computes NET Photosynthetically Active Radiation `PAR` [mol/m²/s].
"""
@inline function compute_PAR(photo::LUEPhotosynthesis{NF}, swdown::NF) where NF
    PAR = NF(0.5) * swdown * (NF(1.0) - photo.α_leaf) * photo.cq
    return PAR
end

"""
    $SIGNATURES
Computes absorbed PAR limited by the fraction of PAR assimilated at ecosystem level `APAR` [mol/m²/s],
Eq. 62, PALADYN (Willeit 2016).
"""
@inline function compute_APAR(photo::LUEPhotosynthesis{NF}, swdown::NF, LAI::NF) where NF
    PAR = compute_PAR(photo, swdown)
    APAR = photo.α_a * PAR * (NF(1.0) - exp(-photo.k_ext*LAI)) 
    return APAR
end

"""
    $SIGNATURES
Computes intercellular CO2 partial pressure [Pa],
Eq. 67, PALADYN (Willeit 2016).
"""
@inline function compute_pres_i(photo::LUEPhotosynthesis, λc, pres_a) 
    pres_i = λc * pres_a
    return pres_i
end

"""
    $SIGNATURES

Computes the temperature stress factor `t_stress` based on the air temperature.
"""

@inline function compute_t_stress(photo::LUEPhotosynthesis{NF}, T_air::NF) where NF
    # TODO check physical meaning of these parameters
    k1 = NF(2.0) * log(NF(1.0)/NF(0.99)-NF(1.0)) / (photo.t_CO2_low - photo.t_photos_low)
    k2 = NF(0.5) * (photo.t_CO2_low + photo.t_photos_low)
    k3 = log(NF(0.99)/NF(0.01)) / (photo.t_CO2_high - photo.t_photos_high)

    # Compute t_stress
    if photo.t_CO2_low < T_air < photo.t_CO2_high
        low = NF(1.0) / (NF(1.0) + exp(k1 * (k2 - T_air)))
        high = NF(1.0) - NF(0.01) * exp(k3 * (T_air - photo.t_photos_high))
        t_stress = low * high
    else
        t_stress = zero(NF)
    end

    return t_stress
end

"""
    $SIGNATURES
Computes factor for light-limited assimilation `c_1` and factor for RuBisCO-limited assimilation `c_2`,
Eqs. C4+C5, PALADYN (Willeit 2016).
"""
@inline function compute_c1_c2(photo::LUEPhotosynthesis{NF}, T_air::NF, Γ_star::NF, Kc::NF, Ko::NF, pres_i::NF, pres_O2::NF) where NF
    t_stress = compute_t_stress(photo, T_air)
    # TODO check factor 2 missing in PALADYN paper
    # TODO add units
    c_1 = photo.α_C3 * t_stress * photo.C_mass * (pres_i - Γ_star) / (pres_i + NF(2.0) * Γ_star)
    c_2 = (pres_i - Γ_star) / (pres_i + Kc * (NF(1.0) + pres_O2 / Ko))
    return c_1, c_2
end

"""
    $SIGNATURES
Computes the maximum rate of net photosynthesis `Vc_max` [gC/m²/s],
following the coordination hypothesis (acclimation), see Harrison 2021 Box 2.
Note: this is not the same formula in PALADYN paper, this implementaion is taken from the code
"""
@inline function compute_Vc_max(photo::LUEPhotosynthesis{NF}, c_1::NF, PAR::NF, Kc::NF, Ko::NF, Γ_star::NF, pres_i::NF, pres_O2::NF) where NF
    Vc_max = c_1 * PAR * (pres_i + Kc * (NF(1.0) + pres_O2 / Ko)) / (pres_i - Γ_star)
    return Vc_max
end

"""
    $SIGNATURES
Computes the PAR-limited and the rubisco-activity-limited photosynthesis rates `JE` and `JC` [gC/m²/s],
Eqn 3+5, Haxeltine & Prentice 1996.
"""
@inline function compute_JE_JC(photo::LUEPhotosynthesis{NF}, c_1::NF, c_2::NF, APAR::NF, Vc_max::NF) where NF
    JE = c_1 * APAR
    JC = c_2 * Vc_max
    return JE, JC
end

"""
    $SIGNATURES
Computes the leaf respiration rate `Rd` [gC/m²/s],
Eqn 10, Haxeltine & Prentice 1996 and Eq. 10 PALADYN (Willeit 2016).
"""
@inline function compute_Rd(photo::LUEPhotosynthesis, Vc_max, β)
    Rd = photo.α_C3 * Vc_max * β
    return Rd
end

"""
    $SIGNATURES
Computes the gross photosynthesis rate `Ag` [gC/m²/s],
Eqn 2, Haxeltine & Prentice 1996
"""
@inline function compute_Ag(photo::LUEPhotosynthesis{NF}, c_1::NF, c_2::NF, APAR::NF, Vc_max::NF, β::NF) where NF
    # Compute JE and Jc, PAR-limited and rubisco-activity-limited photosynthesis rates 
    JE, JC = compute_JE_JC(photo, c_1, c_2, APAR, Vc_max)

    # TODO photosyntheis downregulation ignored for now
    Ag = (JE + JC - sqrt((JE + JC)^2 - NF(4) * photo.θ_r * JE * JC)) / (NF(2) * photo.θ_r) * β
    return Ag
end

"""
    $SIGNATURES

Computes Gross Primary Production `GPP`in [kgC/m²/day] and leaf respiration `Rd` in [gC/m²/day]
"""
function compute_photosynthesis(
    photo::LUEPhotosynthesis{NF},
    T_air::NF, swdown::NF, pres::NF, co2::NF, LAI::NF, λc::NF, β::NF,
) where NF
    # Compute partial CO2 and O2 pressures in [Pa]
    pres_O2 = partial_pressure_O2(pres)
    pres_a = partial_pressure_CO2(pres, co2)

    # TODO check this condition
    if swdown > zero(NF) && T_air > NF(-3.0)
        # Compute kinetic parameters 
        # TODO check physical meaning of these parameters,  Appendix C in PALADYN paper
        # TODO add units
        τ, Kc, Ko = compute_kinetic_parameters(photo, T_air)

        # Compute Γ_star
        # TODO add unit
        Γ_star = compute_Γ_star(photo, τ, pres_O2)

        # TODO check for bioclimatic limit ignored for now
        if LAI > zero(NF)

            # Compute absorbed PAR [mol/m²/s]
            APAR = compute_APAR(photo, swdown, LAI)

            # Compute pres_i, intercellular CO2 partial pressure [Pa]
            pres_i = compute_pres_i(photo, λc, pres_a)
            
            # Compute c1 and c2 parameters for C3 photosynthesis
            # TODO add units
            c_1, c_2 = compute_c1_c2(photo, T_air, Γ_star, Kc, Ko, pres_i, pres_O2)

            # Compute Vc_max, maximum rate of carboxylation [gC/m²/s]
            Vc_max = compute_Vc_max(photo, c_1, APAR, Kc, Ko, Γ_star, pres_i, pres_O2)
            
            # Compute leaf respiration rate [gC/m²/s]
            Rd = compute_Rd(photo, Vc_max, β)

            # Compute Ag, the gross photosynthesis rate [gC/m²/s]
            Ag = compute_Ag(photo, c_1, c_2, APAR, Vc_max, β)

            # Compute An, the net photosynthesis rate [gC/m²/s]
            An = Ag - Rd
        else
            # No leaves, no photosynthesis 
            An = zero(NF)
            # No leaves, no leaf respiration
            Rd = zero(NF)
        end
    else
        # No light, no photosynthesis
        An = zero(NF)
        # TODO Rd = 0 here?
        Rd = zero(NF)
    end

    return Rd, An
end

"""
    $SIGNATURES

Compute the Gross Primary Production [kgC/m²/s].
"""
function compute_GPP(::LUEPhotosynthesis{NF}, An::NF) where {NF}
    # Compute And, the total daytime net photosynthesis `And` [gC/m²/day],
    # Eq 19, Haxeltine & Prentice 1996 + Eq. 65, PALADYN (Willeit 2016).
    # And = (An + (NF(1.0) - day_length / NF(24.0)) * Rd) * seconds_per_day(NF)

    # Compute GPP [kgC/m²/s]
    GPP = An * NF(1.e-3)
    return GPP
end

function compute_auxiliary!(state, model, photo::LUEPhotosynthesis)
    grid = get_grid(model)
    atmos = get_atmosphere(model)
    launch!(grid, :xy, compute_photosynthesis_kernel!, state, photo, atmos)
end

@kernel function compute_photosynthesis_kernel!(
    state,
    photo::LUEPhotosynthesis{NF},
    atmos::AbstractAtmosphere
) where NF
    # TODO checks for positive/negative values in the original PALADYN code ignored for now
    i, j = @index(Global, NTuple)

    # Get inputs
    @inbounds let
        T_air = air_temperature(i, j, state, atmos),
        pres = air_pressure(i, j, state, atmos),
        swdown = shortwave_in(i, j, state, atmos),
        # day_length = daytime_length(i, j, state, atmos),
        co2 = state.CO2[i, j], # no method for this currently...
        β = state.SMLF[i, j],
        LAI = state.LAI[i, j],
        λc = state.λc[i, j];

        # Compute Rd, leaf respiration rate in [gC/m²/s],
        # An, daily net photosynthesis [gC/m²/s]
        Rd, An = compute_photosynthesis(photo, T_air, swdown, pres, co2, LAI, λc, β)

        # Compute GPP, Gross Primary Production in [kgC/m²/s]
        GPP = compute_GPP(photo, An)
    
        # Store results
        state.GPP[i, j] = GPP
        state.Rd[i, j] = Rd
        state.An[i, j] = An
    end
end
