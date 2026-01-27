"""
Return the current value of the `surface_net_radiation` at the given indices.
"""
@propagate_inbounds surface_net_radiation(i, j, grid, state, ::AbstractRadiativeFluxes) = state.surface_net_radiation[i, j]

"""
    $TYPEDEF

Represents the simplest scheme for the radiative budget where outgoing shortwave and longwave
radiation are given as input variables. Net radiation is diagnosed by summing all radiative fluxes:

```math
R_{\\text{net}} = S_{\\uparrow} - S_{\\downarrow} + L_{\\uparrow} - L_{\\downarrow}
```
"""
struct PrescribedRadiativeFluxes <: AbstractRadiativeFluxes end

# Process methods

variables(::PrescribedRadiativeFluxes) = (
    input(:surface_shortwave_up, XY(), units=u"W/m^2", desc="Outgoing (upwelling) shortwave radiation"),
    input(:surface_longwave_up, XY(), units=u"W/m^2", desc="Outgoing (upwelling) longwave radiation"),
    auxiliary(:surface_net_radiation, XY(), units=u"W/m^2", desc="Net outgoing (positive up) radiation"),
)

function compute_auxiliary!(state, model, rad::PrescribedRadiativeFluxes)
    (; grid, atmosphere) = model
    launch!(grid, XY, compute_net_radiation_kernel!, state, rad, atmosphere)
end

"""
    $TYPEDEF

Computes outgoing shortwave and longwave radiation according to separately specified
schemes for the albedo, skin temperature, and atmospheric inputs.
"""
struct DiagnosedRadiativeFluxes <: AbstractRadiativeFluxes end

"""
    compute_surface_net_radiation(::AbstractRadiativeFluxes, surface_shortwave_down, surface_shortwave_up, surface_longwave_down, surface_longwave_up)

Compute the net radiation budget given incoming and outgoing shortwave and longwave radiation.
"""
@inline function compute_surface_net_radiation(::AbstractRadiativeFluxes, surface_shortwave_down, surface_shortwave_up, surface_longwave_down, surface_longwave_up)
    # Sum up radiation fluxes; note that by convention fluxes are positive upward
    surface_net_radiation = surface_shortwave_up - surface_shortwave_down + surface_longwave_up - surface_longwave_down
    return surface_net_radiation
end

"""
    compute_shortwave_out(::DiagnosedRadiativeFluxes, surface_shortwave_down, α)

Compute outgoing shortwave radiation from the incoming `surface_shortwave_down` and albedo `α`.
"""
@inline function compute_shortwave_out(::DiagnosedRadiativeFluxes, surface_shortwave_down, α)
    surface_shortwave_up = α * surface_shortwave_down
    return surface_shortwave_up
end

"""
    compute_longwave_out(::DiagnosedRadiativeFluxes, constants::PhysicalConstants, surface_longwave_down, Ts, ϵ)

Compute outgoing longwave radiation from incoming `surface_longwave_down`, surface temperature `Ts`, and emissivity `ϵ`.
"""
@inline function compute_longwave_out(::DiagnosedRadiativeFluxes, constants::PhysicalConstants, surface_longwave_down, Ts, ϵ)
    T = celsius_to_kelvin(constants, Ts)
    # outgoing LW radiation is the sum of the radiant emittance and the residual incoming radiation
    surface_longwave_up = stefan_boltzmann(constants, T, ϵ) + (1 - ϵ) * surface_longwave_down
    return surface_longwave_up
end

# Process methods

variables(::DiagnosedRadiativeFluxes) = (
    auxiliary(:surface_shortwave_up, XY(), units=u"W/m^2", desc="Outgoing (upwelling) shortwave radiation"),
    auxiliary(:surface_longwave_up, XY(), units=u"W/m^2", desc="Outgoing (upwelling) longwave radiation"),
    auxiliary(:surface_net_radiation, XY(), units=u"W/m^2", desc="Net radiation budget"),
)

function compute_auxiliary!(state, model, rad::DiagnosedRadiativeFluxes)
    (; grid, surface_energy_balance, atmosphere, constants) = model
    launch!(grid, XY, compute_radiative_fluxes_kernel!, state,
            rad, atmosphere, surface_energy_balance.skin_temperature, surface_energy_balance.albedo, constants)
end

# Kernels

@kernel function compute_radiative_fluxes_kernel!(
    state, grid,
    rad::AbstractRadiativeFluxes,
    atmos::AbstractAtmosphere,
    skinT::AbstractSkinTemperature,
    abd::AbstractAlbedo,
    consts::PhysicalConstants
)
    i, j = @index(Global, NTuple)

    # get inputs
    surface_shortwave_up = state.surface_shortwave_up[i, j]
    surface_longwave_up = state.surface_longwave_up[i, j]
    surface_shortwave_down = shortwave_in(i, j, grid, state, atmos)
    surface_longwave_down = longwave_in(i, j, grid, state, atmos)
    Tsurf = skin_temperature(i, j, grid, state, skinT)
    α = albedo(i, j, grid, state, abd)
    ϵ = emissivity(i, j, grid, state, abd)

    # compute outputs
    state.surface_shortwave_up[i, j, 1] = surface_shortwave_up = compute_shortwave_out(rad, surface_shortwave_down, α)
    state.surface_longwave_up[i, j, 1] = surface_longwave_up = compute_longwave_out(rad, consts, surface_longwave_down, Tsurf, ϵ)
    state.surface_net_radiation[i, j, 1] = compute_surface_net_radiation(rad, surface_shortwave_down, surface_shortwave_up, surface_longwave_down, surface_longwave_up)
end

@kernel function compute_net_radiation_kernel!(state, grid, rad::AbstractRadiativeFluxes, atmos::AbstractAtmosphere)
    i, j = @index(Global, NTuple)
    # get inputs
    surface_shortwave_up = state.surface_shortwave_up[i, j]
    surface_longwave_up = state.surface_longwave_up[i, j]
    surface_shortwave_down = shortwave_in(i, j, grid, state, atmos)
    surface_longwave_down = longwave_in(i, j, grid, state, atmos)
    
    # compute net radiation
    state.surface_net_radiation[i, j, 1] = compute_surface_net_radiation(rad, surface_shortwave_down, surface_shortwave_up, surface_longwave_down, surface_longwave_up)
end
