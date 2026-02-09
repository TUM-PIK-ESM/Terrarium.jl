"""
    $TYPEDEF

Canopy interception and storage implementation following PALADYN (Willeit 2016) considering only liquid water (no snow).

Properties:
$FIELDS
"""
@kwdef struct PALADYNCanopyInterception{NF} <: AbstractCanopyInterception{NF}
    "Canopy water interception factor for tree PFTs"
    α_int::NF = 0.2 

    # TODO: Duplicated
    "Extinction coefficient for radiation through vegetation [-]"
    k_ext::NF = 0.5

    "Canopy interception capacity parameter, Verseghy 1991 [kg/m²]"
    w_can_max::NF = 0.2

    "Canopy water removal timescale [s]"
    τ_w::NF = 86400.0
end

PALADYNCanopyInterception(::Type{NF}; kwargs...) where {NF} = PALADYNCanopyInterception{NF}(; kwargs...)

variables(::PALADYNCanopyInterception) = (
    prognostic(:canopy_water, XY(); desc="Canopy liquid water", units=u"kg/m^2"), 
    auxiliary(:canopy_water_interception, XY(); desc="Canopy rain interception rate", units=u"m/s"), 
    auxiliary(:canopy_water_removal, XY(); desc="Canopy water removal rate", units=u"m/s"),
    auxiliary(:saturation_canopy_water, XY(); desc="Fraction of the canopy saturated with water"),
    auxiliary(:precip_ground, XY(); desc="Rainfall rate reaching the ground", units=u"m/s"),
    input(:leaf_area_index, XY(); desc="Leaf Area Index", units=u"m^2/m^2"), 
    input(:SAI, XY(); desc="Stem Area Index", units=u"m^2/m^2"),
)

@propagate_inbounds canopy_water(i, j, grid, fields, ::PALADYNCanopyInterception) = fields.canopy_water[i, j]

@propagate_inbounds saturation_canopy_water(i, j, grid, fields, ::PALADYNCanopyInterception) = fields.saturation_canopy_water[i, j]

@propagate_inbounds ground_precipitation(i, j, grid, fields, ::PALADYNCanopyInterception) = fields.precip_ground[i, j]

"""
    $TYPEDSIGNATURES

Compute `I_can`, the canopy rain interception, following Eq. 42, PALADYN (Willeit 2016).
"""
@inline function compute_canopy_interception(canopy_interception::PALADYNCanopyInterception{NF}, precip, LAI, SAI) where NF   
    I_can = canopy_interception.α_int * precip * (one(NF) - exp(-canopy_interception.k_ext * (LAI + SAI))) 
    return I_can
end

"""
    $TYPEDSIGNATURES

Compute the canopy saturation fraction as `w_can / w_can_max`.
"""
@inline function compute_canopy_saturation_fraction(canopy_interception::PALADYNCanopyInterception{NF}, w_can, LAI, SAI) where NF
    # Compute the wet canopy fraction
    w_can_max = canopy_interception.w_can_max * (LAI + SAI)
    f_can = w_can_max > 0 ? w_can / w_can_max : zero(NF)
    return f_can
end

"""
    $TYPEDSIGNATURES

Compute the canopy water removal rate as `w_can / ρw / τw`.
"""
@inline function compute_canopy_water_removal(
   canopy_interception::PALADYNCanopyInterception{NF},
   constants::PhysicalConstants{NF},
   w_can
) where {NF}
    # Canopy water storage tendency: interception - evaporation - removal
    R_can = max(w_can, zero(NF)) / constants.ρw / canopy_interception.τ_w
    return R_can
end

"""
    $TYPEDSIGNATURES

Compute the `w_can` tendency and removal rate following Eq. 41, PALADYN (Willeit 2016).
"""
@inline function compute_w_can_tendency(
    ::PALADYNCanopyInterception{NF},
    I_can, E_can, R_can
) where NF
    # Canopy water storage tendency: interception - evaporation - removal
    w_can_tend = I_can - E_can - R_can
    return w_can_tend
end

"""
    $SIGNATURES

Compute `precip_ground`, the rate of rain reaching the ground, following a modified version
of Eq. 44, PALADYN (Willeit 2016). Instead of subtracting the tendency, we just directly subtract
interception and add the removal rate `R_can`.
"""
@inline function compute_precip_ground(
    ::PALADYNCanopyInterception{NF},
    precip, I_can, R_can
) where NF
    # Compute the precipitation reaching the ground:
    # precip - interception + removal
    precip_ground = precip - I_can + R_can
    return precip_ground
end

# Process methods

function compute_auxiliary!(
    state, grid,
    canopy_interception::PALADYNCanopyInterception,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    out = auxiliary_fields(state, canopy_interception)
    fields = get_fields(state, canopy_interception, atmos; except = out)
    launch!(grid, XY, compute_auxiliary_kernel!, out, fields, canopy_interception, atmos, constants)
end

function compute_tendencies!(
    state, grid,
    canopy_interception::PALADYNCanopyInterception,
    evtr::AbstractEvapotranspiration,
    args...
)
    out = tendency_fields(state, canopy_interception)
    fields = get_fields(state, canopy_interception, evtr; except = out)
    launch!(grid, XY, compute_tendencies_kernel!, out, fields, canopy_interception, evtr)
end

# Kernel functions

@propagate_inbounds function compute_canopy_auxiliary!(
    out, i, j, grid, fields,
    canopy_interception::PALADYNCanopyInterception{NF},
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants,
    args...
) where NF
    # Get inputs 
    precip = rainfall(i, j, grid, fields, atmos)
    LAI = fields.leaf_area_index[i, j]
    SAI = fields.SAI[i, j]
    w_can = fields.canopy_water[i, j]

    # Compute canopy saturation faction
    f_can = compute_canopy_saturation_fraction(canopy_interception, w_can, LAI, SAI)

    # Compute canopy rain interception
    I_can = compute_canopy_interception(canopy_interception, precip, LAI, SAI)

    # Compute canopy water removal
    R_can = compute_canopy_water_removal(canopy_interception, constants, w_can)

    # Compute precipitation reaching the ground
    precip_ground = compute_precip_ground(canopy_interception, precip, I_can, R_can)

    # Store results
    out.canopy_water_interception[i, j, 1] = I_can
    out.canopy_water_removal[i, j, 1] = R_can
    out.saturation_canopy_water[i, j, 1] = f_can
    out.precip_ground[i, j, 1] = precip_ground
    return out
end

@propagate_inbounds function compute_canopy_water_tendency!(
    tendencies, i, j, grid, fields,
    canopy_interception::PALADYNCanopyInterception,
    evtr::AbstractEvapotranspiration,
    args...
)
    # Get inputs
    E_can = fields.evaporation_canopy[i, j]
    I_can = fields.canopy_water_interception[i, j]
    R_can = fields.canopy_water_removal[i, j]

    # Compute canopy water tendency
    tendencies.canopy_water[i, j, 1] = compute_w_can_tendency(canopy_interception, I_can, E_can, R_can)
    return tendencies
end

# Kernels

@kernel inbounds=true function compute_auxiliary_kernel!(out, grid, fields, canopy_interception::AbstractCanopyInterception, args...)
    i, j = @index(Global, NTuple)
    compute_canopy_auxiliary!(out, i, j, grid, fields, canopy_interception, args...)
end

@kernel inbounds=true function compute_tendencies_kernel!(out, grid, fields, canopy_interception::AbstractCanopyInterception, args...)
    i, j = @index(Global, NTuple)
    compute_canopy_water_tendency!(out, i, j, grid, fields, canopy_interception, args...)
end
