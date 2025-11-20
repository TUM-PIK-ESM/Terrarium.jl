"""
    $TYPEDEF

Canopy hydrology implementation following PALADYN (Willeit 2016) considering only canopy interception and ignoring 
canopy evaporation for now.

Authors: Maha Badri 

Properties:
$TYPEDFIELDS
"""
@kwdef struct PALADYNCanopyHydrology{NF} <: AbstractCanopyHydrology
    
    "Canopy water interception factor"
    α_int::NF = 0.2 

    "Extinction coefficient for radiation"
    k_ext::NF = 0.5

    "Canopy water removal timescale"
    τ_w::NF = 86400.0 # [s]
end

PALADYNCanopyHydrology(::Type{NF}; kwargs...) where {NF} = PALADYNCanopyHydrology{NF}(; kwargs...)

variables(::PALADYNCanopyHydrology) = (
    prognostic(:w_can, XY()), # Canopy liquid water [kg/m²]
    auxiliary(:I_can, XY()), # Canopy rain interception [kg/m²/s]
    input(:Precip, XY()), # Precipitation rate [kg/m²/s]
    input(:LAI, XY()), # Leaf Area Index [m²/m²]
    input(:SAI, XY()) # Stem Area Index [m²/m²]
)


"""
    $SIGNATURES

Computes `I_can`, the canopy rain interception, following Eq. 42, PALADYN (Willeit 2016).
"""
# TODO missing SAI, Precip reahcing the ground
@inline function compute_I_can(canopy_hydrology::PALADYNCanopyHydrology{NF}, Precip_ground, LAI, SAI) where NF   
    I_can = canopy_hydrology.α_int * Precip_ground * (1 - exp(-canopy_hydrology.k_ext*(LAI + SAI))) 
    return I_can
end

"""
    $SIGNATURES
Computes the `w_can` tendency following Eq. 41, PALADYN (Willeit 2016).
Modified to ignore the evaporation term for now.
"""

@inline function compute_w_can_tend(canopy_hydrology::PALADYNCanopyHydrology{NF}, Precip_ground, LAI, SAI, w_can) where NF
    # Compute I_can
    I_can = compute_I_can(canopy_hydrology, Precip_ground, LAI, SAI)

    # Compute w_can tendency
    w_can_tendency = I_can - w_can / canopy_hydrology.τ_w

    return w_can_tendency
end

function compute_auxiliary!(state, model, canopy_hydrology::PALADYNCanopyHydrology)
    grid = get_grid(model)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, canopy_hydrology)
end

@kernel function compute_auxiliary_kernel!(state, canopy_hydrology::PALADYNCanopyHydrology{NF}) where NF
    i, j = @index(Global, NTuple)

    # Compute canopy rain interception
    # Assume all precipitation reaches the ground for now
    Precip_ground = state.Precip[i, j]
    LAI = state.LAI[i, j]
    SAI = state.SAI[i, j]
    state.I_can[i, j] = compute_I_can(canopy_hydrology, Precip_ground, LAI, SAI)
end

function compute_tendencies!(state, model, canopy_hydrology::PALADYNCanopyHydrology)
    grid = get_grid(model)
    launch!(grid, :xy, compute_tendencies_kernel!, state, canopy_hydrology)
end

@kernel function compute_tendencies_kernel!(state, canopy_hydrology::PALADYNCanopyHydrology{NF}) where NF  
    i, j = @index(Global, NTuple)

    # Get inputs
    # Assume all precipitation reaches the ground for now
    Precip_ground = state.Precip[i, j]
    LAI = state.LAI[i, j]
    SAI = state.SAI[i, j]
    w_can = state.w_can[i, j]

    # Compute the canopy water tendency
    w_can_tendency = compute_w_can_tend(canopy_hydrology, Precip_ground, LAI, SAI, w_can)
    
    # Store result
    state.tendencies.w_can[i, j] = w_can_tendency
end