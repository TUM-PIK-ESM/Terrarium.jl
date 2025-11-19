# Note: maybe change the name later, if the PALADYN autotrophic respiration approach has a more specific name
"""
    $TYPEDEF

Autotrophic respiration implementation from PALADYN (Willeit 2016).

Authors: Maha Badri and Matteo Willeit

Properties:
$TYPEDFIELDS
"""
@kwdef struct PALADYNAutotrophicRespiration{NF} <: AbstractAutotrophicRespiration
    # TODO check physical meaning of this parameter + add unit
    "Sapwood parameter"
    cn_sapwood::NF = 330.0

    # TODO check physical meaning of this parameter + add unit
    "Root parameter"
    cn_root::NF = 29.0

    "Ratio of total to respiring stem carbon, Cox 2001, PFT specific [-]"
    aws::NF = 10.0 # Value for Needleleaf tree PFT
end

variables(::PALADYNAutotrophicRespiration) = (
    auxiliary(:Rd, XY()), # Daily leaf respiration [gC/m²/day]
    auxiliary(:phen, XY()), # Phenology factor [-]
    auxiliary(:GPP, XY()), # Gross Primary Production [kgC/m²/day]
    auxiliary(:Ra, XY()), # Autotrophic respiration [kgC/m²/day]
)

"""
    $SIGNATURES

Computes temperature factors `f_temp_air` and `f_temp_soil` for autotrophic respiration.
"""

@inline function compute_f_temp(
        autoresp::PALADYNAutotrophicRespiration{NF},
        T_air::NF,
    ) where {NF}
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
@inline function compute_resp10(autoresp::PALADYNAutotrophicRespiration{NF}) where {NF}
    # TODO check physical meaning of this variable + add unit
    # TODO add resp10 implementation
    # For now, placeholder as a constant value
    resp10 = NF(0.066)

    return resp10
end

"""
$SIGNATURES

Computes maintenance respiration `Rm` in [kgC/m²/day].
"""
@inline function compute_Rm(
        autoresp::PALADYNAutotrophicRespiration{NF},
        vegcarbon_dynamics::PALADYNCarbonDynamics{NF},
        T_air,
        Rd,
        phen,
        C_veg,
    ) where {NF}

    # Compute f_temp for autotrophic respiration
    f_temp_air, f_temp_soil = compute_f_temp(autoresp, T_air)

    # Compute resp10
    resp10 = compute_resp10(autoresp)

    # Compute leaf respiration
    R_leaf = Rd / NF(1000.0) # convert from gC/m²/day to kgC/m²/day

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

Computes growth respiration `Rg` in [kgC/m²/day].
"""
@inline function compute_Rg(autoresp::PALADYNAutotrophicRespiration{NF}, GPP, Rm) where {NF}
    Rg = NF(0.25) * (GPP - Rm)
    return Rg
end

"""
$SIGNATURES

Computes autotrophic respiration `Ra` as the sum of maintenance respiration `Rm` and growth respiration `Rg` in [kgC/m²/day].
"""
@inline function compute_Ra(autoresp::PALADYNAutotrophicRespiration, vegcarbon_dynamics::PALADYNCarbonDynamics, T_air, Rd, phen, C_veg, GPP)
    # Compute Rm, maintenance respiration
    Rm = compute_Rm(autoresp, vegcarbon_dynamics, T_air, Rd, phen, C_veg)

    # Compute Rg, growth respiration
    Rg = compute_Rg(autoresp, GPP, Rm)

    # Compute Ra, autotrophic respiration
    Ra = Rm + Rg
    return Ra
end

"""
$SIGNATURES

Computes Net Primary Productivity `NPP` as the difference between Gross Primary Production `GPP` and autotrophic respiration `Ra`
in [kgC/m²/day].
"""
@inline function compute_NPP(autoresp::PALADYNAutotrophicRespiration, GPP, Ra)
    NPP = GPP - Ra
    return NPP
end


function compute_auxiliary!(state, model, autoresp::PALADYNAutotrophicRespiration)
    grid = get_grid(model)
    bcs = get_boundary_conditions(model)
    return launch!(grid, :xy, compute_auxiliary_kernel!, state, autoresp, get_carbon_dynamics(model), bcs.top)
end

@kernel function compute_auxiliary_kernel!(
        state,
        autoresp::PALADYNAutotrophicRespiration{NF},
        vegcarbon_dynamics::PALADYNCarbonDynamics{NF},
        ::PrescribedAtmosphere
    ) where {NF}
    i, j = @index(Global, NTuple)

    # Get inputs
    T_air = state.T_air[i, j]
    Rd = state.Rd[i, j]
    phen = state.phen[i, j]
    C_veg = state.C_veg[i, j]
    GPP = state.GPP[i, j]

    # Compute autotrophic respiration Ra
    Ra = compute_Ra(autoresp, vegcarbon_dynamics, T_air, Rd, phen, C_veg, GPP)

    # Compute NPP
    NPP = compute_NPP(autoresp, GPP, Ra)

    # Store results
    state.Ra[i, j] = Ra
    state.NPP[i, j] = NPP
end
