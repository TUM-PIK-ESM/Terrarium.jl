
# Note: photosynthesis implementation following PALADYN for C3 PFTs

@kwdef struct LUEPhotosynthesis{NF} <: AbstractPhotosynthesis
    # TODO check physical meaning of this parameter 
    "Value of τ at 25°C"
    τ25::NF = 2600.0 

    # TODO check physical meaning of this parameter 
    "Value of kc at 25°C"
    Kc25::NF = 30.0

    # TODO check physical meaning of this parameter 
    "Value of ko at 25°C"
    Ko25::NF = 3.0e4

    # TODO check physical meaning of this parameter
    "q10 for temperature-sensitive parameter τ"
    q10_τ::NF = 0.57

    # TODO check physical meaning of this parameter
    "q10 for temperature-sensitive parameter kc"
    q10_Kc::NF = 2.1

    # TODO check physical meaning of this parameter
    "q10 for temperature-sensitive parameter ko"
    q10_Ko::NF = 1.2

    "Leaf albedo in PAR range"
    α_leaf::NF = 0.17

    "Conversion factor for solar radiation at 550 nm from J/m² to mol/m² [mol/J]"
    cq::NF = 4.6e-6

    "Extinction coefficient for radiation through vegetation"
    k_ext::NF = 0.5

    "Fraction of PAR assimilated at ecosystem level, relative to leaf level"
    αa::NF = 0.5

    #TODO check physical meaning of this parameter
    "Parameter, PFT specific"
    t_CO2_high::NF = 42 # Value for Needleleaf tree PFT 

    # TODO check physical meaning of this parameter
    "Parameter, PFT specific"
    t_CO2_low::NF = -4.0 # Value for Needleleaf tree PFT 

    # TODO check physical meaning of this parameter
    "Parameter, PFT specific"
    t_photos_high::NF = 30 # Value for Needleleaf tree PFT 

    # TODO check physical meaning of this parameter
    "Parameter, PFT specific"
    t_photos_low::NF = 15 # Value for Needleleaf tree PFT 

    "Intrinsic quantum efficiency of CO2 uptake in C3 plants"
    α_C3::NF = 0.08 

    "Atomic mass of carbon"
    C_mass::NF = 12.0

    # TODO check physical meaning of this parameter
    "Shape parameter"
    θr::NF = 0.7 

end

variables(photo::LUEPhotosynthesis) = (
    # TODO for now define atmospheric inputs/forcings here, move later
    auxiliary(:T_air, XY()), # Surface air temperature in Kelvin [K]
    auxiliary(:q_air, XY()), # Surface air specific humidity [kg/kg]
    auxiliary(:p, XY()), # Surface pressure [Pa]
    auxiliary(:swdown, XY()), # Downwelling shortwave radiation at the surface [W/m²]
    auxiliary(:co2, XY()), # Atmospheric CO2 concentration [ppm]
    auxiliary(:GPP, XY()), # Gross Primary Production [kgC/m²/day]
)

# TODO for now define functions that compute derived variables from atm. inputs/forcings here, move later
@inline function convert_T_to_Celsius(idx, state, model::AbstractVegetationModel, photo::LUEPhotosynthesis) 
    i, j = idx

    # Get model constants
    constants = get_constants(model)

    # Convert T_air to °C
    T_air_C = state.T_air[i, j] - constants.T0

    return T_air_C
end

# TODO for now define functions that compute derived variables from atm. inputs/forcings here, move later
@inline function compute_p_O2(idx, state, model::AbstractVegetationModel, photo::LUEPhotosynthesis{NF}) where NF
    i, j = idx

    # Get surface pressure
    p = state.p[i, j]

    # Compute O2 partial pressure [Pa]
    p_O2 = NF(0.209) * p

    return p_O2
end

# TODO for now define functions that compute derived variables from atm. inputs/forcings here, move later
@inline function compute_pa(idx, state, model::AbstractVegetationModel, photo::LUEPhotosynthesis{NF}) where NF
    i, j = idx
    
    # Get surface pressure
    p = state.p[i, j]
    
    # Compute atmospheric CO2 partial pressure [Pa]
    pa = co2 * NF(1e-6) * p

    return pa
end

@inline function compute_f_temp(idx, state, model::AbstractVegetationModel, photo::LUEPhotosynthesis{NF}) where NF
    i, j = idx

    # Convert surface air temperature to Celsius
    T_air_C = convert_T_to_Celsius(idx, state, model, photo)

    # TODO check physical meaning of these parameters
    k1 = NF(2.0) * log(NF(1.0)/NF(0.99)-NF(1.0)) / (photo.t_CO2_low - photo.t_photos_low)
    k2 = NF(0.5) * (photo.t_CO2_low + photo.t_photos_low)
    k3 = log(NF(0.99)/NF(0.01)) / (photo.t_CO2_high - photo.t_photos_high)

    # Compute f_temp, a PFT-specific temperature inhibition function
    if T_air_C < photo.t_CO2_high && T_air_C > photo.t_CO2_low
        low = NF(1.0) / (NF(1.0) + exp(k1 * (k2 - T_air_C)))
        high = NF(1.0) - NF(0.01) * exp(k3 * (T_air_C - photo.t_photos_high))
        f_temp = low * high
    else
        f_temp = NF(0.0)
    end

    return f_temp
    
end

@inline function compute_kinetic_parameters(idx, state, model::AbstractVegetationModel, photo::LUEPhotosynthesis{NF}) where NF
     i, j = idx

    # Convert surface air temperature to Celsius
    T_air_C = convert_T_to_Celsius(idx, state, model, photo)

    # TODO check meaning of these parameters, Appendix C in PALADYN paper
    τ = photo.τ25 * photo.q10_τ^((T_air_C - NF(25.0)) * NF(0.1))
    Kc = photo.Kc25 * photo.q10_Kc^((T_air_C - NF(25.0)) * NF(0.1))
    Ko = photo.Ko25 * photo.q10_Ko^((T_air_C - NF(25.0)) * NF(0.1))
    # TODO is Γ_star a kinetic parameter?
    p_O2 = compute_p_O2(idx, state, model, photo)
    Γ_star = p_O2 / (NF(2.0) * τ)
    
    return τ, Kc, Ko, Γ_star
end

@inline function compute_auxiliary!(idx, state, model::AbstractVegetationModel, photo::LUEPhotosynthesis{NF}) where NF
    # TODO checks for positive/negative values in the original PALADYN code ignored for now

    i, j = idx

    # Get atmospheric inputs/forcings and compute derived variables
    swdown = state.swdown[i, j] 
    T_air_C = convert_T_to_Celsius(idx, state, model, photo)
    p_O2 = compute_p_O2(idx, state, model, photo)
    pa = compute_pa(idx, state, model, photo)

    # TODO add daylength/sec_day implementaion
    # For now, placeholders as constant values
    daylength = NF(24.0)
    sec_day = NF(8.765813e4)

    # TODO check this condition
    if (daylength > NF(0.0)) && (T_air_C > NF(-3.0))

        # Compute kinetic parameters 
        # TODO check physical meaning of these parameters,  Appendix C in PALADYN paper
        τ, Kc, Ko, Γ_star = compute_kinetic_parameters(idx, state, model, photo)

        # Compute NET photosynthetically active radiation [mol/m²/day]
        par = NF(0.5) * swdown * sec_day * (NF(1.0) - photo.α_leaf) * photo.cq

        # TODO check for bioclimatic limit ignored for now
        if state.LAI[i, j] > NF(0.0)
            # Compute absorbed PAR limited by the fraction of PAR assimilated at ecosystem level, the leaf scattering
            apar = photo.αa * (NF(1.0) - exp(-photo.k_ext*state.LAI[i, j])) * par
            
            # Compute intercellular CO2 partial pressure
            pi = state.λc[i, j] * pa

            # Compute c1 and c2 parameters for C3 photosynthesis
            # TODO check factor 2 missing in PALADYN paper
            c_1 = photo.α_C3 * f_temp * photo.C_mass * (pi - Γ_star) / (pi + NF(2.0) * Γ_star) 
            c_2 = (pi - Γ_star) / (pi + Kc*(NF(1.0) + p_O2/Ko))

            # Compute the maximum daily rate of net photosynthesis [gC/m²/day]
            # Following the coordination hypothesis (acclimation), see Harrison 2021 Box 2
            # Note: this is not the same formula in PALADYN paper, this implementaion is taken from the code
            Vc_max = c_1 * apar * (pi + Kc * (NF(1.0) + p_O2 / Ko)) / (pi - Γ_star)

            # Compute the PAR-limited photosynthesis rate [molC/m²/h]
            # Eqn 3, Haxeltine & Prentice 1996
            JE = c_1 * apar / daylength

            # Compute the rubisco-activity-limited photosynthesis rate [molC/m²/h]
            # Eqn 5, Haxeltine & Prentice 1996
            JC = c_2 * Vc_max / NF(24.0)

            # TODO add implementaion for the soil-moisture limiting factor
            # For now, set it to 1.0, no soil moisture limitation
            β = NF(1.0)

            # Compute the daily gross photosynthesis [gC/m²/day]
            # Eqn 2, Haxeltine & Prentice 1996
            # TODO photosyntheis downregulation ignored for now
            Ag = (JE + JC - sqrt((JE + JC)^2 - NF(4.0) * θr * JE * JC)) / (NF(2.0) * θr) * daylength * β

            # Compute the daily leaf respiration [gC/m2/day]
            # Eqn 10, Haxeltine & Prentice 1996
            Rd = α_C3 * Vc_max * β

            # Compute the daily net photosynthesis [gC/m²/day]
            An = Ag - Rd

            # Compute total daytime net photosynthesis [gC/m²/day]
            # Eqn 19, Haxeltine & Prentice 1996
            And = An + (NF(1.0) - daylength / NF(24.0)) * Rd

            # Compute daily GPP [kgC/m²/day]
            state.GPP[i, j] = And * NF(1.e-3)
        else
            # No photosynthesis 
            state.GPP[i, j] = NF(0.0)
        end
    else
        # No light
        state.GPP[i, j] = NF(0.0)
    end

end
