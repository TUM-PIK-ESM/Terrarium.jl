"""
Default implementation of [net_radiation](@ref) for all `AbstractRadiativeFluxes` that
simply returns the current value of the `RadNet` field.
"""
@inline net_radiation(idx, state, ::AbstractRadiativeFluxes) = state.RadNet[idx...]

struct PrescribedRadiativeFluxes <: AbstractRadiativeFluxes end

variables(::PrescribedRadiativeFluxes) = (
    input(:SwOut, XY(), units=u"W/m^2", desc="Outoing (upwelling) shortwave radiation"),
    input(:LwOut, XY(), units=u"W/m^2", desc="Outoing (upwelling) longwave radiation"),
    auxiliary(:RadNet, XY(), units=u"W/m^2", desc="Net radiation budget"),
)

function compute_auxiliary!(state, model, rad::PrescribedRadiativeFluxes)
    grid = get_grid(model)
    atmos = get_atmosphere(model)
    launch!(grid, compute_net_radiation!, state, rad, atmos)
end

struct DiagnosedRadiativeFluxes <: AbstractRadiativeFluxes end

variables(::DiagnosedRadiativeFluxes) = (
    auxiliary(:SwOut, XY(), units=u"W/m^2", desc="Outoing (upwelling) shortwave radiation"),
    auxiliary(:LwOut, XY(), units=u"W/m^2", desc="Outoing (upwelling) longwave radiation"),
    auxiliary(:RadNet, XY(), units=u"W/m^2", desc="Net radiation budget"),
)

function compute_auxiliary!(state, model, rad::DiagnosedRadiativeFluxes)
    grid = get_grid(model)
    atmos = get_atmosphere(model)
    constants = get_constants(model)
    launch!(grid, compute_radiative_fluxes!, state, rad, atmos, constants)
end

# Kernels

@kernel function compute_radiative_fluxes!(
    state, grid,
    rad::AbstractRadiativeFluxes,
    atmos::AbstractAtmosphere,
    skinT::AbstractSkinTemperature,
    albedo::AbstractAlbedo,
    constants::PhysicalConstants
)
    i, j = @index(Global, NTuple)
    idx = (i, j)

    # get inputs
    SwOut = state.SwOut[i, j]
    LwOut = state.LwOut[i, j]
    SwIn = shortwave_in(idx, state, atmos)
    LwIn = longwave_in(idx, state, atmos)
    Tsurf = skin_temperature(idx, state, skinT)
    α = albedo(idx, state, albedo)
    ϵ = emissivity(idx, state, albedo)

    # compute outputs
    state.SwOut[i, j, 1] = SwOut = shortwave_out(rad, SwIn, α)
    state.LwOut[i, j, 1] = LwOut = longwave_out(rad, constants, LwIn, Tsurf, ϵ)
    state.RadNet[i, j, 1] = net_radiation(rad, SwIn, SwOut, LwIn, LwOut)
end

@kernel function compute_net_radiation!(state, rad::AbstractRadiativeFluxes, atmos::AbstractAtmosphere)
    i, j = @index(Global, NTuple)
    # get inputs
    SwOut = state.SwOut[i, j]
    LwOut = state.LwOut[i, j]
    SwIn = shortwave_in(idx, state, atmos)
    LwIn = longwave_in(idx, state, atmos)
    
    # compute net radiation
    state.RadNet[i, j, 1] = net_radiation(rad, SwIn, SwOut, LwIn, LwOut)
end

# Kernel functions

"""
    net_radiation(::AbstractRadiativeFluxes, SwIn, SwOut, LwIn, LwOut)

Compute the net radiation budget given incoming and outgoing shortwave and longwave radiation.
"""
@inline function net_radiation(::AbstractRadiativeFluxes, SwIn, SwOut, LwIn, LwOut)
    # Sum up radiation fluxes; note that by convention fluxes are positive upward
    RadNet = SwOut - SwIn + LwOut - LwIn
    return RadNet
end

"""
    shortwave_out(::DiagnosedRadiativeFluxes, SwIn, α)

Compute outgoing shortwave radiation from the incoming `SwIn` and albedo `α`.
"""
@inline function shortwave_out(::DiagnosedRadiativeFluxes, SwIn, α)
    SwOut = α * SwIn
    return SwOut
end

"""
    longwave_out(::DiagnosedRadiativeFluxes, constants::PhysicalConstants, LwIn, Ts, ϵ)

Compute outgoing longwave radiation from incoming `LwIn`, surface temperature `Ts`, and emissivity `ϵ`.
"""
@inline function longwave_out(::DiagnosedRadiativeFluxes, constants::PhysicalConstants, LwIn, Ts, ϵ)
    T = celsius_to_kelvin(constants, Ts)
    # outgoing LW radiation is the sum of the radiant emittance and the residual incoming radiation
    LwOut = stefan_boltzmann(constants, T, ϵ) + (1 - ϵ) * LwIn
    return LwOut
end
