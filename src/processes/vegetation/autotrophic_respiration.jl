# Note: maybe change the name later, if the PALADYN autotrophic respiration approach has a more specific name
"""
    $TYPEDEF

Autotrophic respiration implementation from PALADYN (Willeit 2016).

Authors: Maha Badri and Matteo Willeit

Properties:
$TYPEDFIELDS
"""
@kwdef struct PALADYNAutotrophicRespiration{NF} <: AbstractAutotrophicRespiration
    # TODO check physical meaning of this parameter
    "Sapwood parameter"
    cn_sapwood::NF = 330.0

    # TODO check physical meaning of this parameter
    "Root parameter"
    cn_root::NF = 29.0

    "Ratio of total to respiring stem carbon, Cox 2001, PFT specific"
    aws::NF = 10.0 # Value for Needleleaf tree PFT
end

variables(autoresp::PALADYNAutotrophicRespiration) = (
    auxiliary(:T_air, XY()), # Surface air temperature in Celcius [°C]
    auxiliary(:GPP, XY()), # Gross Primary Production [kgC/m²/day]
    auxiliary(:Rd, XY()), # Respiration at 10°C [kgC/m²/day]
    auxiliary(:C_veg, XY()), # Vegetation carbon pool [kgC/m²]
    auxiliary(:phen, XY()), # Phenology factor
    auxiliary(:Ra, XY()), # Autotrophic respiration [kgC/m²/day]
    auxiliary(:NPP, XY()), # Net Primary Production [kgC/m²/day]
)


"""
    $SIGNATURES

Computes temperature factors `f_temp_air` and `f_temp_soil` for autotrophic respiration.
"""

@inline function compute_f_temp(
    autoresp::PALADYNAutotrophicRespiration{NF},
    T_air::NF,
) where NF
    # Compute f_temp_soil
    # TODO add f_temp_soil implementaion (depends on soil temperature)
    # For now, placeholder as a constant value
    f_temp_soil = zero(NF)

    # Compute f_temp_air
    # TODO: These hardcoded constants need to be moved either into the model struct as
    # parameters or into the PhysicalConstants struct
    f_temp_air = exp(NF(308.56) * (NF(1.0) / NF(56.02) - NF(1.0) / (NF(46.02) + T_air)))

    return f_temp_air, f_temp_soil
end

"""
$SIGNATURES

Computes `resp10` 
"""
@inline function compute_resp10(autoresp::PALADYNAutotrophicRespiration{NF}) where NF
    # TODO check physical meaning of this variable
    # TODO add resp10 implementation
    # For now, placeholder as a constant value
    resp10 = NF(0.066)

    return resp10
end

"""
$SIGNATURES

Computes maintenance respiration `Rm`.
"""
@inline function compute_Rm(
    autoresp::PALADYNAutotrophicRespiration{NF}, 
    vegcarbon_dynamics::PALADYNCarbonDynamics{NF}, 
    T_air,
    Rd, 
    phen,
    C_veg,
) where NF

    # Compute f_temp for autotrophic respiration
    f_temp_air, f_temp_soil = compute_f_temp(autoresp, T_air)

    # Compute resp10
    resp10 = compute_resp10(autoresp)
    
    # Compute leaf respiration
    R_leaf = Rd/NF(1000.0) # Convert from gC/m²/day to kgC/m²/day

    # Compute stem respiration
    R_stem = resp10 * f_temp_air * (vegcarbon_dynamics.awl * ((NF(2.0) / vegcarbon_dynamics.SLA) + vegcarbon_dynamics.awl)) / 
                                   (C_veg * autoresp.aws * autoresp.cn_sapwood) 

    # Compute root respiration
    R_root = resp10 * f_temp_soil * phen * (NF(2.0) / vegcarbon_dynamics.SLA) / 
                                                       (vegcarbon_dynamics.SLA * C_veg * autoresp.cn_root)

    # Compute maintenance respiration Rm
    Rm = R_leaf + R_stem + R_root 

    return Rm
end

"""
$SIGNATURES

Computes growth respiration `Rg`.
"""
@inline function compute_Rg(autoresp::PALADYNAutotrophicRespiration{NF}, GPP, Rm) where NF
    Rg = NF(0.25) * (GPP - Rm)
    return Rg
end

"""
$SIGNATURES

Computes autotrophic respiration `Ra` as the sum of maintenance respiration `Rm` and growth respiration `Rg`.
"""
@inline function compute_Ra(autoresp::PALADYNAutotrophicRespiration, Rm, Rg) 
    Ra = Rm + Rg
    return Ra
end

"""
$SIGNATURES

Computes Net Primary Productivity `NPP` as the difference between Gross Primary Production `GPP` and autotrophic respiration `Ra`.
"""
@inline function compute_NPP(autoresp::PALADYNAutotrophicRespiration, GPP, Ra) 
    NPP = GPP - Ra
    return NPP
end


function compute_auxiliary!(state, model, autoresp::PALADYNAutotrophicRespiration)
    grid = get_grid(model)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, autoresp, get_carbon_dynamics(model))
end

@kernel function compute_auxiliary_kernel!(
    state,
    autoresp::PALADYNAutotrophicRespiration{NF},
    vegcarbon_dynamics::PALADYNCarbonDynamics{NF},
) where NF
    i, j = @index(Global, NTuple)
    
    # Compute maintenance respiration
    Rm = compute_Rm(autoresp, vegcarbon_dynamics, state.T_air[i, j], state.Rd[i, j], state.phen[i, j], state.C_veg[i, j])

    # Compute growth respiration Rg
    Rg = compute_Rg(autoresp, state.GPP[i, j], Rm)

    # Compute autotrophic respiration Ra
    state.Ra[i, j] = compute_Ra(autoresp, Rm, Rg)

    # Compute NPP
    state.NPP[i, j] = compute_NPP(autoresp, state.GPP[i, j], state.Ra[i, j])
end

   