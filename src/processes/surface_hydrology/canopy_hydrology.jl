"""
    $TYPEDEF

Canopy interception and storage implementation following PALADYN (Willeit 2016) considering only liquid water (no snow).

Properties:
$FIELDS
"""
@kwdef struct PALADYNCanopyHydrology{NF} <: AbstractCanopyHydrology
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

PALADYNCanopyHydrology(::Type{NF}; kwargs...) where {NF} = PALADYNCanopyHydrology{NF}(; kwargs...)

variables(::PALADYNCanopyHydrology) = (
    prognostic(:canopy_water, XY(); desc="Canopy liquid water", units=u"kg/m^2"), 
    auxiliary(:canopy_water_interception, XY(); desc="Canopy rain interception rate", units=u"m/s"), 
    auxiliary(:canopy_water_removal, XY(); desc="Canopy water removal rate", units=u"m/s"),
    auxiliary(:saturation_canopy_water, XY(); desc="Fraction of the canopy saturated with water"),
    auxiliary(:precip_ground, XY(); desc="Rainfall rate reaching the ground", units=u"m/s"),
    input(:LAI, XY(); desc="Leaf Area Index", units=u"m^2/m^2"), 
    input(:SAI, XY(); desc="Stem Area Index", units=u"m^2/m^2"),
)

@propagate_inbounds canopy_water(i, j, state, grid, ::PALADYNCanopyHydrology) = state.canopy_water[i, j]

@propagate_inbounds saturation_canopy_water(i, j, state, grid, ::PALADYNCanopyHydrology) = state.saturation_canopy_water[i, j]

@propagate_inbounds ground_precipitation(i, j, state, grid, ::PALADYNCanopyHydrology) = state.precip_ground[i, j]

"""
    $TYPEDSIGNATURES

Compute `I_can`, the canopy rain interception, following Eq. 42, PALADYN (Willeit 2016).
"""
@inline function compute_canopy_interception(canopy_hydrology::PALADYNCanopyHydrology{NF}, precip, LAI, SAI) where NF   
    I_can = canopy_hydrology.α_int * precip * (one(NF) - exp(-canopy_hydrology.k_ext*(LAI + SAI))) 
    return I_can
end

"""
    $TYPEDSIGNATURES

Compute the canopy saturation fraction as `w_can / w_can_max`.
"""
@inline function compute_canopy_saturation_fraction(canopy_hydrology::PALADYNCanopyHydrology{NF}, w_can, LAI, SAI) where NF
    # Compute the wet canopy fraction
    w_can_max = canopy_hydrology.w_can_max * (LAI + SAI)
    f_can = w_can_max > 0 ? w_can / w_can_max : zero(NF)
    return f_can
end

"""
    $TYPEDSIGNATURES

Compute the canopy water removal rate as `w_can / ρw / τw`.
"""
@inline function compute_canopy_water_removal(
   canopy_hydrology::PALADYNCanopyHydrology{NF},
   constants::PhysicalConstants{NF},
   w_can
) where {NF}
    # Canopy water storage tendency: interception - evaporation - removal
    R_can = max(w_can, zero(NF)) / constants.ρw / canopy_hydrology.τ_w
    return R_can
end

"""
    $TYPEDSIGNATURES

Compute the `w_can` tendency and removal rate following Eq. 41, PALADYN (Willeit 2016).
"""
@inline function compute_w_can_tend(
    ::PALADYNCanopyHydrology{NF},
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
    ::PALADYNCanopyHydrology{NF},
    precip, I_can, R_can
) where NF
    # Compute the precipitation reaching the ground:
    # precip - interception + removal
    precip_ground = precip - I_can + R_can
    return precip_ground
end

# Process methods

function compute_auxiliary!(state, model, canopy_hydrology::PALADYNCanopyHydrology)
    grid = get_grid(model)
    atmos = get_atmosphere(model)
    constants = get_constants(model)
    launch!(state, grid, :xy, compute_auxiliary_kernel!, canopy_hydrology, atmos, constants)
end

function compute_tendencies!(state, model, canopy_hydrology::PALADYNCanopyHydrology)
    grid = get_grid(model)
    launch!(state, grid, :xy, compute_tendencies_kernel!, canopy_hydrology)
end

# Kernels

@kernel inbounds=true function compute_auxiliary_kernel!(
    state, 
    canopy_hydrology::PALADYNCanopyHydrology{NF},
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
) where NF
    
    i, j = @index(Global, NTuple)

    # Get inputs 
    precip = rainfall(i, j, state, grid, atmos)
    LAI = state.LAI[i, j]
    SAI = state.SAI[i, j]
    w_can = state.canopy_water[i, j]

    # Compute canopy saturation faction
    f_can = compute_canopy_saturation_fraction(canopy_hydrology, w_can, LAI, SAI)

    # Compute canopy rain interception
    I_can = compute_canopy_interception(canopy_hydrology, precip, LAI, SAI)

    # Compute canopy water removal
    R_can = compute_canopy_removal(canopy_hydrology, constants, w_can)

    # Compute precipitation reaching the ground
    precip_ground = compute_precip_ground(canopy_hydrology, precip, I_can, R_can)

    # Store results
    state.canopy_water_interception[i, j] = I_can
    state.canopy_water_removal[i, j] = R_can
    state.saturation_canopy_water[i, j] = f_can
    state.precip_ground[i, j] = precip_ground
end

@kernel function compute_tendencies_kernel!(
    state, 
    canopy_hydrology::PALADYNCanopyHydrology{NF},
    canopy_ET::PALADYNCanopyEvapotranspiration{NF}
) where NF  
    i, j = @index(Global, NTuple)

    # Get inputs
    @inbounds let I_can = state.canopy_water_interception[i, j],
                  E_can = state.evaporation_canopy[i, j],
                  R_can = state.canopy_water_removal;

        # Compute canopy water tendency
        state.tendencies.w_can[i, j] = compute_w_can_tend(canopy_hydrology, I_can, E_can, R_can)
    end
end
