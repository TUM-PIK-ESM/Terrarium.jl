"""
    $TYPEDEF

Canopy interception and storage implementation following PALADYN (Willeit 2016) considering only liquid water (no snow).

Properties:
$TYPEDFIELDS
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
    auxiliary(:saturation_water_canopy, XY(); desc="Fraction of the canopy saturated with water"),
    auxiliary(:precip_ground, XY(); desc="Rainfall rate reaching the ground", units=u"m/s"),
    input(:LAI, XY(); desc="Leaf Area Index", units=u"m^2/m^2"), 
    input(:SAI, XY(); desc="Stem Area Index", units=u"m^2/m^2"),
)

@inline @propagate_inbounds saturation_water_canopy(i, j, state, grid, ::PALADYNCanopyHydrology) = state.saturation_water_canopy[i, j]

function compute_auxiliary!(state, model, canopy_hydrology::PALADYNCanopyHydrology)
    grid = get_grid(model)
    atmos = get_atmosphere(model)
    constants = get_constants(model)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, canopy_hydrology, atmos, constants)
end

function compute_tendencies!(state, model, canopy_hydrology::PALADYNCanopyHydrology)
    grid = get_grid(model)
    launch!(grid, :xy, compute_tendencies_kernel!, state, canopy_hydrology)
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
    state.saturation_water_canopy[i, j] = f_can
end

@kernel function compute_tendencies_kernel!(
    state, 
    canopy_hydrology::PALADYNCanopyHydrology{NF},
    canopy_ET::CanopyEvapotranspiration{NF}
) where NF  
    i, j = @index(Global, NTuple)

    # Get inputs
    @inbounds let w_can = state.canopy_water[i, j],
                  I_can = state.interception_water_canopy[i, j],
                  E_can = state.evaporation_canopy[i, j];

        # Compute canopy water tendency
        state.tendencies.w_can[i, j] = w_can_tend = compute_w_can_tend(canopy_hydrology, w_can, I_can, E_can)

        # Compute precipitation reaching the ground
        state.precip_ground[i, j] = compute_precip_ground(canopy_hydrology, precip, E_can, w_can_tend)
    end
end


# Kernel functions

"""
    $SIGNATURES

Compute `I_can`, the canopy rain interception, following Eq. 42, PALADYN (Willeit 2016).
"""
@inline function compute_canopy_interception(canopy_hydrology::PALADYNCanopyHydrology{NF}, precip, LAI, SAI) where NF   
    I_can = canopy_hydrology.α_int * precip * (one(NF) - exp(-canopy_hydrology.k_ext*(LAI + SAI))) 
    return I_can
end

"""
    $SIGNATURES

Compute the canopy saturation fraction as `w_can / w_can_max`.
"""
@inline function compute_canopy_saturation_fraction(canopy_hydrology::PALADYNCanopyHydrology, w_can, LAI, SAI)
    # Compute the wet canopy fraction
    w_can_max = canopy_hydrology.can_max_w * (LAI + SAI)
    f_can = w_can / w_can_max
    return f_can
end

"""
    $SIGNATURES
Compute the `w_can` tendency following Eq. 41, PALADYN (Willeit 2016).
"""
@inline function compute_w_can_tend(canopy_hydrology::PALADYNCanopyHydrology{NF}, I_can, E_can, R_can) where NF
    # Canopy water storage tendency: interception - evaporation - runoff
    w_can_tend = I_can - E_can - R_can
    return w_can_tend
end

"""
    $SIGNATURES

Compute `precip_ground`, the rate of rain reaching the ground, following Eq. 44, PALADYN (Willeit 2016).
"""
@inline function compute_precip_ground(
    canopy_hydrology::PALADYNCanopyHydrology{NF},
    constants::PhysicalConstants{NF},
    precip, E_can, w_can_tend
) where NF
    # Compute the precipitation reaching the ground
    precip_ground = precip - E_can - w_can_tend / constants.ρw
    return precip_ground
end
