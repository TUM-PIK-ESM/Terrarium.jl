# Surface energy balance processes

"""
Base type for surface energy balance schemes.
"""
abstract type AbstractSurfaceEnergyBalance{NF} <: AbstractProcess end

"""
    compute_surface_energy_fluxes!(state, grid, ::AbstractSurfaceEnergyBalance, args...)

Compute the surface energy fluxes and skin temperature from the current `state` and `grid`.
The required `args` are implementation dependent.
"""
function compute_surface_energy_fluxes! end

"""
Base type for skin temperature and ground heat flux schemes.
"""
abstract type AbstractSkinTemperature{NF} <: AbstractProcess end

"""
    skin_temperature(i, j, grid, state, ::AbstractSkinTemperature)

Return the current skin temperature at the given indices.
"""
@propagate_inbounds skin_temperature(i, j, grid, state, ::AbstractSkinTemperature) = state.skin_temperature[i, j]

"""
    ground_heat_flux(i, j, grid, state, ::AbstractSkinTemperature)

Return the current ground heat flux at the given indices.
"""
@propagate_inbounds ground_heat_flux(i, j, grid, state, ::AbstractSkinTemperature) = state.ground_heat_flux[i, j]

"""
    $TYPEDEF

Base type for radiative flux parameterizations.
"""
abstract type AbstractRadiativeFluxes{NF} <: AbstractProcess end

"""
    shortwave_up(i, j, grid, state, ::AbstractRadiativeFluxes)

Return the current outgoing (upwelling) shortwave radiation at the surface.
"""
@propagate_inbounds shortwave_up(i, j, grid, state, ::AbstractRadiativeFluxes) = state.shortwave_up[i, j]

"""
    longwave_up(i, j, grid, state, ::AbstractRadiativeFluxes)

Return the current outgoing (upwelling) longwave radiation at the given indices `i, j`.
"""
@propagate_inbounds longwave_up(i, j, grid, state, ::AbstractRadiativeFluxes) = state.shortwave_up[i, j]

"""
    surface_net_radiation(i, j, grid, state, ::AbstractRadiativeFluxes)

Return the current surface net radiation at the given indices `i, j`.
"""
@propagate_inbounds surface_net_radiation(i, j, grid, state, ::AbstractRadiativeFluxes) = state.surface_net_radiation[i, j]

"""
Base type for turbulent (latent and sensible) heat flux parameterizations.
"""
abstract type AbstractTurbulentFluxes{NF} <: AbstractProcess end

"""
    sensible_heat_flux(i, j, grid, state, ::AbstractTurbulentFluxes)

Return the current sensible heat flux at the given indices.
"""
@propagate_inbounds sensible_heat_flux(i, j, grid, state, ::AbstractRadiativeFluxes) = state.sensible_heat_flux[i, j]

"""
    latent_heat_flux(i, j, grid, state, ::AbstractTurbulentFluxes)

Return the current latent heat flux at the given indices.
"""
@propagate_inbounds latent_heat_flux(i, j, grid, state, ::AbstractRadiativeFluxes) = state.latent_heat_flux[i, j]

# Surface energy balance parameterizations

"""
Base type for surface albedo and emissivity parameterizations.
"""
abstract type AbstractAlbedo{NF} end

"""
    albedo(i, j, grid, state, ::AbstractAlbedo)

Return the current albedo at the given indices.
"""
@propagate_inbounds albedo(i, j, grid, state, ::AbstractAlbedo) = state.albedo[i, j]

"""
    emissivity(i, j, grid, state, ::AbstractAlbedo)

Return the current emissivity at the given indices.
"""
@propagate_inbounds emissivity(i, j, grid, state, ::AbstractAlbedo) = state.emissivity[i, j]
