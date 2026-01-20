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
    auxiliary(:interception_water_canopy, XY(); desc="Canopy rain interception", units=u"m/s"), 
    auxiliary(:saturation_canopy_water, XY(); desc="Fraction of the canopy saturated with water"),
    auxiliary(:precip_ground, XY(); desc="Rainfall rate reaching the ground", units=u"m/s"),
    input(:LAI, XY(); desc="Leaf Area Index", units=u"m^2/m^2"), 
    input(:SAI, XY(); desc="Stem Area Index", units=u"m^2/m^2"),
)

@inline @propagate_inbounds canopy_water(i, j, state, grid, ::PALADYNCanopyHydrology) = state.canopy_water[i, j]

@inline @propagate_inbounds saturation_canopy_water(i, j, state, grid, ::PALADYNCanopyHydrology) = state.saturation_canopy_water[i, j]

@inline @propagate_inbounds ground_precipitation(i, j, state, grid, ::PALADYNCanopyHydrology) = state.precip_ground[i, j]

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

@kernel function compute_auxiliary_kernel!(
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

    # Store results
    state.interception_water_canopy[i, j] = I_can
    state.saturation_canopy_water[i, j] = f_can
end

@kernel function compute_tendencies_kernel!(
    state, 
    canopy_hydrology::PALADYNCanopyHydrology{NF},
    canopy_ET::PALADYNCanopyEvapotranspiration{NF}
) where NF  
    i, j = @index(Global, NTuple)

    # Get inputs
    @inbounds let w_can = state.canopy_water[i, j],
                  I_can = state.interception_water_canopy[i, j],
                  E_can = state.evaporation_canopy[i, j];

        # Compute canopy water tendency
        w_can_tend, R_can = compute_w_can_tend(canopy_hydrology, constants, w_can, I_can, E_can)
        state.tendencies.w_can[i, j] = w_can_tend

        # Compute precipitation reaching the ground
        state.precip_ground[i, j] = compute_precip_ground(canopy_hydrology, constants, precip, I_can, R_can)
    end
end


# Kernel functions

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

Compute the `w_can` tendency and removal rate following Eq. 41, PALADYN (Willeit 2016).
"""
@inline function compute_w_can_tend(
    canopy_hydrology::PALADYNCanopyHydrology{NF},
    constants::PhysicalConstants{NF},
    w_can, I_can, E_can
) where NF
    # Canopy water storage tendency: interception - evaporation - removal
    R_can = w_can / constants.ρw / canopy_hydrology.τ_w
    w_can_tend = I_can - E_can - R_can
    return w_can_tend, R_can
end

"""
    $SIGNATURES

Compute `precip_ground`, the rate of rain reaching the ground, following a modified version
of Eq. 44, PALADYN (Willeit 2016). Instead of subtracting the tendency, we just directly subtract
interception and add the removal rate `R_can`.
"""
@inline function compute_precip_ground(
    canopy_hydrology::PALADYNCanopyHydrology{NF},
    constants::PhysicalConstants{NF},
    precip, I_can, R_can
) where NF
    # Compute the precipitation reaching the ground:
    # precip - interception + removal
    precip_ground = precip - I_can + R_can
    return precip_ground
end
