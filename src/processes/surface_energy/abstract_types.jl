# Surface energy balance types

"""
Base type for surface albedo and emissivity process implementations.
"""
abstract type AbstractAlbedo <: AbstractProcess end

"""
    albedo(i, j, state, grid, ::AbstractAlbedo)

Return the current albedo at the given indices.
"""
function albedo end

"""
    emissivity(i, j, state, grid, ::AbstractAlbedo)

Return the current emissivity at the given indices.
"""
function emissivity end

"""
Base type for radiation budget schemes.
"""
abstract type AbstractRadiativeFluxes <: AbstractProcess end

"""
    surface_net_radiation(i, j, state, grid, ::AbstractRadiativeFluxes)

Return the current net radiation at the given indices.
"""
function surface_net_radiation end

"""
Base type for turbulent (latent and sensible) heat fluxes at the surface.
"""
abstract type AbstractTurbulentFluxes <: AbstractProcess end

"""
    sensible_heat_flux(i, j, state, grid, ::AbstractTurbulentFluxes)

Return the current sensible heat flux at the given indices.
"""
function sensible_heat_flux end

"""
    latent_heat_flux(i, j, state, grid, ::AbstractTurbulentFluxes)

Return the current latent heat flux at the given indices.
"""
function latent_heat_flux end

"""
Base type for skin temperature and ground heat flux schemes.
"""
abstract type AbstractSkinTemperature <: AbstractProcess end

"""
    skin_temperature(i, j, state, grid, ::AbstractSkinTemperature)

Return the current skin temperature at the given indices.
"""
function skin_temperature end

"""
Base type for surface energy balance schemes.
"""
abstract type AbstractSurfaceEnergyBalance <: AbstractProcess end
