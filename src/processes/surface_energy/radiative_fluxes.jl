"""
Default implementation of [net_incoming_radiation](@ref) for all `AbstractRadiativeFluxes` that
simply returns the current value of the `net_incoming_radiation` field.
"""
@inline net_incoming_radiation(idx, state, ::AbstractRadiativeFluxes) = state.net_incoming_radiation[idx...]

"""
    $TYPEDEF

Represents the simplest scheme for the radiative budget where outgoing shortwave and longwave
radiation are given as input variables. Net radiation is diagnosed by summing all radiative fluxes:

```math
R_{\\text{net}} = S_{\\uparrow} - S_{\\downarrow} + L_{\\uparrow} - L_{\\downarrow}
```
"""
struct PrescribedRadiativeFluxes <: AbstractRadiativeFluxes end

variables(::PrescribedRadiativeFluxes) = (
    input(:surface_shortwave_up, XY(), units=u"W/m^2", desc="Outoing (upwelling) shortwave radiation"),
    input(:surface_longwave_out, XY(), units=u"W/m^2", desc="Outoing (upwelling) longwave radiation"),
    auxiliary(:net_incoming_radiation, XY(), units=u"W/m^2", desc="Net radiation budget"),
)

function compute_auxiliary!(state, model, rad::PrescribedRadiativeFluxes)
    (; grid, atmosphere) = model
    launch!(state, grid, :xy, compute_net_radiation!, rad, atmosphere)
end

"""
    $TYPEDEF

Computes outgoing shortwave and longwave radiation according to separately specified
schemes for the albedo, skin temperature, and atmospheric inputs.
"""
struct DiagnosedRadiativeFluxes <: AbstractRadiativeFluxes end

variables(::DiagnosedRadiativeFluxes) = (
    auxiliary(:surface_shortwave_up, XY(), units=u"W/m^2", desc="Outoing (upwelling) shortwave radiation"),
    auxiliary(:surface_longwave_out, XY(), units=u"W/m^2", desc="Outoing (upwelling) longwave radiation"),
    auxiliary(:net_incoming_radiation, XY(), units=u"W/m^2", desc="Net radiation budget"),
)

function compute_auxiliary!(state, model, rad::DiagnosedRadiativeFluxes)
    (; grid, surface_energy_balance, atmosphere, constants) = model
    launch!(
        state,
        grid,
        :xy,
        compute_radiative_fluxes!,
        rad,
        atmosphere,
        surface_energy_balance.skin_temperature,
        surface_energy_balance.albedo,
        constants
    )
end

# Kernels

@kernel function compute_radiative_fluxes!(
    state,
    ::AbstractLandGrid,
    rad::AbstractRadiativeFluxes,
    atmos::AbstractAtmosphere,
    skinT::AbstractSkinTemperature,
    albd::AbstractAlbedo,
    consts::PhysicalConstants
)
    i, j = @index(Global, NTuple)
    idx = (i, j)

    # get inputs
    surface_shortwave_up = state.surface_shortwave_up[i, j]
    surface_longwave_out = state.surface_longwave_out[i, j]
    surface_shortwave_down = shortwave_in(idx, state, atmos)
    surface_longwave_down = longwave_in(idx, state, atmos)
    Tsurf = skin_temperature(idx, state, skinT)
    α = albedo(idx, state, albd)
    ϵ = emissivity(idx, state, albd)

    # compute outputs
    state.surface_shortwave_up[i, j, 1] = surface_shortwave_up = shortwave_out(rad, surface_shortwave_down, α)
    state.surface_longwave_out[i, j, 1] = surface_longwave_out = longwave_out(rad, consts, surface_longwave_down, Tsurf, ϵ)
    state.net_incoming_radiation[i, j, 1] = net_incoming_radiation(rad, surface_shortwave_down, surface_shortwave_up, surface_longwave_down, surface_longwave_out)
end

@kernel function compute_net_radiation!(state, ::AbstractLandGrid, rad::AbstractRadiativeFluxes, atmos::AbstractAtmosphere)
    i, j = @index(Global, NTuple)
    idx = (i, j)
    # get inputs
    surface_shortwave_up = state.surface_shortwave_up[i, j]
    surface_longwave_out = state.surface_longwave_out[i, j]
    surface_shortwave_down = shortwave_in(idx, state, atmos)
    surface_longwave_down = longwave_in(idx, state, atmos)
    
    # compute net radiation
    state.net_incoming_radiation[i, j, 1] = net_incoming_radiation(rad, surface_shortwave_down, surface_shortwave_up, surface_longwave_down, surface_longwave_out)
end

# Kernel functions

"""
    net_incoming_radiation(::AbstractRadiativeFluxes, surface_shortwave_down, surface_shortwave_up, surface_longwave_down, surface_longwave_out)

Compute the net radiation budget given incoming and outgoing shortwave and longwave radiation.
"""
@inline function net_incoming_radiation(::AbstractRadiativeFluxes, surface_shortwave_down, surface_shortwave_up, surface_longwave_down, surface_longwave_out)
    # Sum up radiation fluxes; note that by convention fluxes are positive upward
    net_incoming_radiation = surface_shortwave_up - surface_shortwave_down + surface_longwave_out - surface_longwave_down
    return net_incoming_radiation
end

"""
    shortwave_out(::DiagnosedRadiativeFluxes, surface_shortwave_down, α)

Compute outgoing shortwave radiation from the incoming `surface_shortwave_down` and albedo `α`.
"""
@inline function shortwave_out(::DiagnosedRadiativeFluxes, surface_shortwave_down, α)
    surface_shortwave_up = α * surface_shortwave_down
    return surface_shortwave_up
end

"""
    longwave_out(::DiagnosedRadiativeFluxes, constants::PhysicalConstants, surface_longwave_down, Ts, ϵ)

Compute outgoing longwave radiation from incoming `surface_longwave_down`, surface temperature `Ts`, and emissivity `ϵ`.
"""
@inline function longwave_out(::DiagnosedRadiativeFluxes, constants::PhysicalConstants, surface_longwave_down, Ts, ϵ)
    T = celsius_to_kelvin(constants, Ts)
    # outgoing LW radiation is the sum of the radiant emittance and the residual incoming radiation
    surface_longwave_out = stefan_boltzmann(constants, T, ϵ) + (1 - ϵ) * surface_longwave_down
    return surface_longwave_out
end
