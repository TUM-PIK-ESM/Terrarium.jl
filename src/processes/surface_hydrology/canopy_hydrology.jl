"""
    $TYPEDEF

Canopy hydrology implementation following PALADYN (Willeit 2016) considering only liquid water (no snow).

Properties:
$TYPEDFIELDS
"""
@kwdef struct PALADYNCanopyHydrology{NF} <: AbstractCanopyHydrology
    
    "Canopy water interception factor for tree PFTs"
    α_int::NF = 0.2 

    "Extinction coefficient for radiation"
    k_ext::NF = 0.5

    "Canopy interception capacity parameter, Verseghy 1991 [kg/m²]"
    can_max_w::NF = 0.2

    "Canopy water removal timescale [s]"
    τ_w::NF = 86400.0
end

PALADYNCanopyHydrology(::Type{NF}; kwargs...) where {NF} = PALADYNCanopyHydrology{NF}(; kwargs...)

variables(::PALADYNCanopyHydrology) = (
    prognostic(:w_can, XY(); desc="Canopy liquid water", units=u"kg/m^2"), 
    auxiliary(:I_can, XY(); desc="Canopy rain interception", units=u"kg/m^2/s"), 
    auxiliary(:E_can, XY(); desc="Canopy evaporation", units=u"kg/m^2/s"),
    auxiliary(:precip_ground, XY(); desc="Rainfall rate reaching the ground", units=u"kg/m^2/s"),
    input(:LAI, XY(); desc="Leaf Area Index", units=u"m^2/m^2"), 
    input(:SAI, XY(); desc="Stem Area Index", units=u"m^2/m^2")
)

"""
    $SIGNATURES

Computes `I_can`, the canopy rain interception, following Eq. 42, PALADYN (Willeit 2016).
"""
@inline function compute_I_can(canopy_hydrology::PALADYNCanopyHydrology{NF}, precip, LAI, SAI) where NF   
    I_can = canopy_hydrology.α_int * precip * (one(NF) - exp(-canopy_hydrology.k_ext*(LAI + SAI))) 
    return I_can
end

"""
    $SIGNATURES

Computes `E_can`, the canopy evaporation, following Eq. 43, PALADYN (Willeit 2016).
"""
@inline function compute_E_can(
    canopy_hydrology::PALADYNCanopyHydrology{NF}, 
    constants::PhysicalConstants{NF},   
    w_can, LAI, SAI, rₐ, q_sat, q_air
    ) where NF 
    # Compute the wet canopy fraction
    w_can_max = canopy_hydrology.can_max_w * (LAI + SAI) 
    f_can = w_can / w_can_max
    
    # Compute the canopy evaporation
    E_can = (constants.ρₐ / rₐ) * (q_sat - q_air) * f_can
    return E_can
end


"""
    $SIGNATURES
Computes the `w_can` tendency following Eq. 41, PALADYN (Willeit 2016).
"""
@inline function compute_w_can_tend(canopy_hydrology::PALADYNCanopyHydrology{NF}, w_can, I_can, E_can) where NF
    
    w_can_tend = I_can - E_can - w_can / canopy_hydrology.τ_w

    return w_can_tend
end

"""
    $SIGNATURES

Computes `precip_ground`, the rate of rain reaching the ground, following Eq. 44, PALADYN (Willeit 2016).
"""
@inline function compute_precip_ground(canopy_hydrology::PALADYNCanopyHydrology{NF}, precip, w_can, I_can, E_can) where NF   
    # Compute the canopy water tendency
    w_can_tend = compute_w_can_tend(canopy_hydrology, w_can, I_can, E_can)

    # Compute the precipitation reaching the ground
    precip_ground = precip - E_can - w_can_tend
    return precip_ground
end


function compute_auxiliary!(state, model, canopy_hydrology::PALADYNCanopyHydrology)
    grid = get_grid(model)
    atmos = get_atmosphere(model)
    seb = get_surface_energy_balance(model)
    constants = get_constants(model)
    launch!(grid, :xy, compute_auxiliary_kernel!, state, canopy_hydrology, constants, atmos, seb)
end

@kernel function compute_auxiliary_kernel!(
    state, 
    canopy_hydrology::PALADYNCanopyHydrology{NF},
    constants::PhysicalConstants{NF},
    atmos::AbstractAtmosphere,
    seb::AbstractSurfaceEnergyBalance           
) where NF
    
    i, j = @index(Global, NTuple)

    # Get inputs 
    precip = rainfall(i, j, state, atmos)
    q_air = specific_humidity(i, j, state, atmos)
    rₐ = aerodynamic_resistance(i, j, state, seb.turbulent_fluxes.aerodynamic_resistance)
    q_sat = specific_humidity_saturation(i, j, state, seb.skin_temperature, atmos, constants)
    LAI = state.LAI[i, j]
    SAI = state.SAI[i, j]
    w_can = state.w_can[i, j]

    # Compute canopy rain interception
    I_can = compute_I_can(canopy_hydrology, precip, LAI, SAI)
    
    # Compute canopy evaporation
    E_can = compute_E_can(canopy_hydrology, constants, w_can, LAI, SAI, rₐ, q_sat, q_air)

    # Compute precipitation reaching the ground
    precip_ground = compute_precip_ground(canopy_hydrology, precip, w_can, I_can, E_can)

    # Store results
    state.I_can[i, j] = I_can
    state.E_can[i, j] = E_can
    state.precip_ground[i, j] = precip_ground

end

function compute_tendencies!(state, model, canopy_hydrology::PALADYNCanopyHydrology)
    grid = get_grid(model)
    launch!(grid, :xy, compute_tendencies_kernel!, state, canopy_hydrology)
end

@kernel function compute_tendencies_kernel!(
    state, 
    canopy_hydrology::PALADYNCanopyHydrology{NF}
) where NF  
    i, j = @index(Global, NTuple)

    # Get inputs
    w_can = state.w_can[i, j]
    I_can = state.I_can[i, j]
    E_can = state.E_can[i, j]

    # Compute canopy water tendency
    state.tendencies.w_can[i, j] = compute_w_can_tend(canopy_hydrology, w_can, I_can, E_can)
end