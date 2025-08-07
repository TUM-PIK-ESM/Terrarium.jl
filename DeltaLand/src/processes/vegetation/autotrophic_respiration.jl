# Note: maybe change the name later, if the PALADYN autotrophic respiration approach has a more specific name

@kwdef struct PALADYNAutotrophicRespiration{NF} <: AbstractAutotrophicRespiration
    # TODO check physical meaning of this parameter
    "Parameter"
    cn_sapwood::NF = 330

    # TODO check physical meaning of this parameter
    "Parameter"
    cn_root::NF = 29

    "Ratio of total to respiring stem carbon, Cox 2001, PFT specific"
    aws::NF = 10 # Value for Needleleaf tree PFT
end

variables(autoresp::PALADYNAutotrophicRespiration) = (
    auxiliary(:T_air, XY()), # Surface air temperature [K]
    auxiliary(:GPP, XY()), # Gross Primary Production [kgC/m²/day]
    auxiliary(:Rd, XY()), # Respiration at 10°C [kgC/m²/day]
    auxiliary(:C_veg, XY()), # Vegetation carbon pool [kgC/m²]
    auxiliary(:phen, XY()), # Phenology factor
    auxiliary(:Ra, XY()), # Autotrophic respiration [kgC/m²/day]
    auxiliary(:NPP, XY()), # Net Primary Production [kgC/m²/day]
)

@inline function compute_f_temp(idx, state, model::AbstractVegetationModel, autoresp::PALADYNAutotrophicRespiration{NF}) where NF
    i, j = idx

    # Get atmospheric inputs/forcings and compute derived variables
    T_air_C = convert_T_to_Celsius(idx, state, model, model.photosynthesis)

    # Compute f_temp_soil
    # TODO add f_temp_soil implementaion (depends on soil temperature)
    # For now, placeholder as a constant value
    f_temp_soil = NF(0.0) 

    # Compute f_temp_air
    ftemp_air = exp(NF(308.56) * (NF(1.0) / NF(56.02) - NF(1.0) / (NF(46.02) + T_air_C)))

    return f_temp_air, f_temp_soil
end

@inline function compute_auxiliary!(idx, state, model::AbstractVegetationModel, autoresp::PALADYNAutotrophicRespiration{NF}) where NF
    i, j = idx
    
    # Get carbon_dynamics
    carbon_dynamics = get_carbon_dynamics(model)

    # Compute f_temp for autotrophic respiration
    f_temp_air, f_temp_soil = compute_f_temp(idx, state, model, autoresp)

    # TODO add resp10 implementation
    # TODO check physical meaning of this variable
    # For now, placeholder as a constant value
    resp10 = NF(0.066)

    # Compute leaf respiration
    R_leaf = state.Rd[i, j]/NF(1000.0) # Convert from gC/m²/day to kgC/m²/day

    # Compute stem respiration
    R_stem = resp10 * f_temp_air * (carbon_dynamics.awl * ((NF(2.0) / carbon_dynamics.SLA) + carbon_dynamics.awl)) / 
                                   (state.C_veg[i, j] * autoresp.aws * autoresp.cn_sapwood) 

    # Compute root respiration
    R_root = resp10 * f_temp_soil * state.phen[i, j] * (NF(2.0) / carbon_dynamics.SLA) / 
                                                       (carbon_dynamics.SLA * state.C_veg[i, j] * autoresp.cn_root)

    
    # Compute maintenance respiration Rm
    Rm = R_leaf + R_stem + R_root 

    # Compute growth respiration Rg
    Rg = NF(0.25) * (state.GPP[i, j] - Rm)

    # Compute autotrophic respiration Ra
    state.Ra[i, j] = Rm + Rg

    # Compute NPP
    state.NPP[i, j] = state.GPP[i, j] - state.Ra[i, j]
end

   