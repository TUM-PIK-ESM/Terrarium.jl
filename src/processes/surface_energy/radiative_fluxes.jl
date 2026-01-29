# Prescribed radiative fluxes

"""
    $TYPEDEF

Represents the simplest scheme for the radiative budget where outgoing shortwave and longwave
radiation are given as input variables. Net radiation is diagnosed by summing all radiative fluxes:

```math
R_{\\text{net}} = S_{\\uparrow} - S_{\\downarrow} + L_{\\uparrow} - L_{\\downarrow}
```
"""
struct PrescribedRadiativeFluxes{NF} <: AbstractRadiativeFluxes{NF} end

PrescribedRadiativeFluxes(::Type{NF}) where {NF} = PrescribedRadiativeFluxes{NF}()

## Process methods

variables(::PrescribedRadiativeFluxes) = (
    input(:surface_shortwave_up, XY(), units=u"W/m^2", desc="Outgoing (upwelling) shortwave radiation"),
    input(:surface_longwave_up, XY(), units=u"W/m^2", desc="Outgoing (upwelling) longwave radiation"),
    auxiliary(:surface_net_radiation, XY(), units=u"W/m^2", desc="Net outgoing (positive up) radiation"),
)

@inline function compute_auxiliary!(state, grid, rad::PrescribedRadiativeFluxes)
    launch!(grid, XY, compute_radiative_fluxes_kernel!, state, rad)
end

## Kernel functions

@propagate_inbounds function compute_surface_upwelling_radiation(i, j, grid, state, rad::PrescribedRadiativeFluxes, args...)
    surface_shortwave_up = state.surface_shortwave_up[i, j]
    surface_longwave_up = state.surface_longwave_up[i, j]
    return (; surface_shortwave_up, surface_longwave_up)
end

# Diagnosed radiative fluxes

"""
    $TYPEDEF

Computes outgoing shortwave and longwave radiation according to separately specified
schemes for the albedo, skin temperature, and atmospheric inputs.
"""
struct DiagnosedRadiativeFluxes{NF} <: AbstractRadiativeFluxes{NF} end

DiagnosedRadiativeFluxes(::Type{NF}) where {NF} = DiagnosedRadiativeFluxes{NF}()

"""
    compute_shortwave_up(::DiagnosedRadiativeFluxes, surface_shortwave_down, α)

Compute outgoing shortwave radiation from the incoming `surface_shortwave_down` and albedo `α`.
"""
@inline function compute_shortwave_up(::DiagnosedRadiativeFluxes, surface_shortwave_down, α)
    surface_shortwave_up = α * surface_shortwave_down
    return surface_shortwave_up
end

"""
    compute_longwave_up(::DiagnosedRadiativeFluxes, constants::PhysicalConstants, surface_longwave_down, Ts, ϵ)

Compute outgoing longwave radiation from incoming `surface_longwave_down`, surface temperature `Ts`, and emissivity `ϵ`.
"""
@inline function compute_longwave_up(::DiagnosedRadiativeFluxes, constants::PhysicalConstants, surface_longwave_down, Ts, ϵ)
    T = celsius_to_kelvin(constants, Ts)
    # outgoing LW radiation is the sum of the radiant emittance and the residual incoming radiation
    surface_longwave_up = stefan_boltzmann(constants, T, ϵ) + (1 - ϵ) * surface_longwave_down
    return surface_longwave_up
end

## Process methods

variables(::DiagnosedRadiativeFluxes) = (
    auxiliary(:surface_shortwave_up, XY(), units=u"W/m^2", desc="Outgoing (upwelling) shortwave radiation"),
    auxiliary(:surface_longwave_up, XY(), units=u"W/m^2", desc="Outgoing (upwelling) longwave radiation"),
    auxiliary(:surface_net_radiation, XY(), units=u"W/m^2", desc="Net radiation budget"),
)

@inline function compute_auxiliary!(
    state, grid,
    rad::DiagnosedRadiativeFluxes,
    skinT::AbstractSkinTemperature,
    albedo::AbstractAlbedo,
    atmos::AbstractAtmosphere,
    consts::PhysicalConstants
)
    launch!(grid, XY, compute_radiative_fluxes_kernel!, state, rad, skinT, albedo, atmos, consts)
end

## Kernel functions

@propagate_inbounds function compute_surface_upwelling_radiation(
    i, j, grid, state,
    rad::DiagnosedRadiativeFluxes,
    skinT::AbstractSkinTemperature,
    abd::AbstractAlbedo,
    atmos::AbstractAtmosphere,
    consts::PhysicalConstants
)
    # Get inputs
    SW_down = shortwave_down(i, j, grid, state, atmos)
    LW_down = longwave_down(i, j, grid, state, atmos)
    Tsurf = skin_temperature(i, j, grid, state, skinT)
    α = albedo(i, j, grid, state, abd)
    ϵ = emissivity(i, j, grid, state, abd)

    # Compute fluxes
    surface_shortwave_up = compute_shortwave_up(rad, SW_down, α)
    surface_longwave_up = compute_longwave_up(rad, consts, LW_down, Tsurf, ϵ)
    
    # Return fluxes as named tuple
    return (; surface_shortwave_up, surface_longwave_up)
end

# Common implementations

"""
    $TYPEDEF

Compute the net radiation budget given incoming and outgoing shortwave and longwave radiation.
"""
@inline function compute_surface_net_radiation(
    ::AbstractRadiativeFluxes,
    SW_up,
    SW_down,
    LW_up,
    LW_down
)
    # Sum up radiation fluxes; note that by convention fluxes are positive upward
    rad_net = SW_up - SW_down + LW_up - LW_down
    return rad_net
end

## Kernels

@kernel inbounds=true function compute_radiative_fluxes_kernel!(state, grid, rad::AbstractRadiativeFluxes, args...)
    i, j = @index(Global, NTuple)

    # Compute and store outgoing fluxes
    outgoing_fluxes = compute_surface_upwelling_radiation(i, j, grid, state, rad, args...)
    state.surface_shortwave_up[i, j, 1] = outgoing_fluxes.surface_shortwave_up[i, j]
    state.surface_longwave_up[i, j, 1] = outgoing_fluxes.surface_longwave_up[i, j]
    
    # Comptue and store net radiation
    state.surface_net_radiation[i, j, 1] = compute_surface_net_radiation(i, j, grid, state, rad)
end

## Kernel functions

@propagate_inbounds function compute_surface_net_radiation(i, j, grid, state, rad::AbstractRadiativeFluxes)
    # Get inputs
    SW_down = shortwave_down(i, j, grid, state, atmos)
    SW_up = shortwave_up(i, j, grid, state, rad)
    LW_down = longwave_down(i, j, grid, state, atmos)
    LW_up = longwave_up(i, j, grid, state, rad)
    
    # Compute and return net radiation
    surface_net_radiation = compute_surface_net_radiation(rad, SW_up, SW_down, LW_up, LW_down)
    return surface_net_radiation
end
