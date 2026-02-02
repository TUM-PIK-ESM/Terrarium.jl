# Surface energy balance processes

"""
Base type for surface energy balance schemes which couple together the relevant processes
for radiative and turbulent surface energy fluxes.
"""
abstract type AbstractSurfaceEnergyBalance{NF} <: AbstractComponent{NF} end

## Getter methods for SEB types
"""
    get_radiative_fluxes(seb)

Return the `radiative_fluxes` component of the surface energy balance.
"""
@inline get_radiative_fluxes(seb::AbstractSurfaceEnergyBalance) = seb.radiative_fluxes

"""
    get_turbulent_fluxes(seb)

Return the `turbulent_fluxes` component of the surface energy balance.
"""
@inline get_turbulent_fluxes(seb::AbstractSurfaceEnergyBalance) = seb.turbulent_fluxes

"""
    get_skin_temperature(seb)

Return the `skin_temperature` process from the surface energy balance.
"""
@inline get_skin_temperature(seb::AbstractSurfaceEnergyBalance) = seb.skin_temperature

"""
    get_albedo(seb)

Return the `albedo` parameterization associated with the surface energy balance.
"""
@inline get_albedo(seb::AbstractSurfaceEnergyBalance) = seb.albedo

"""
    compute_surface_energy_fluxes!(state, grid, ::AbstractSurfaceEnergyBalance, args...)

Compute the surface energy fluxes and skin temperature from the current `state` and `grid`.
The required `args` are implementation dependent.
"""
function compute_surface_energy_fluxes! end

"""
Base type for skin temperature and ground heat flux schemes.
"""
abstract type AbstractSkinTemperature{NF} <: AbstractProcess{NF} end

"""
    skin_temperature(i, j, grid, fields, ::AbstractSkinTemperature)

Return the current skin temperature at the given indices.
"""
@propagate_inbounds skin_temperature(i, j, grid, fields, ::AbstractSkinTemperature) = fields.skin_temperature[i, j]

"""
    ground_heat_flux(i, j, grid, fields, ::AbstractSkinTemperature)

Return the current ground heat flux at the given indices.
"""
@propagate_inbounds ground_heat_flux(i, j, grid, fields, ::AbstractSkinTemperature) = fields.ground_heat_flux[i, j]

"""
    $TYPEDEF

Base type for radiative flux parameterizations.
"""
abstract type AbstractRadiativeFluxes{NF} <: AbstractProcess{NF} end

"""
    shortwave_up(i, j, grid, fields, ::AbstractRadiativeFluxes)

Return the current outgoing (upwelling) shortwave radiation at the surface.
"""
@propagate_inbounds shortwave_up(i, j, grid, fields, ::AbstractRadiativeFluxes) = fields.surface_shortwave_up[i, j]

"""
    longwave_up(i, j, grid, fields, ::AbstractRadiativeFluxes)

Return the current outgoing (upwelling) longwave radiation at the given indices `i, j`.
"""
@propagate_inbounds longwave_up(i, j, grid, fields, ::AbstractRadiativeFluxes) = fields.surface_longwave_up[i, j]

"""
    surface_net_radiation(i, j, grid, fields, ::AbstractRadiativeFluxes)

Return the current surface net radiation at the given indices `i, j`.
"""
@propagate_inbounds surface_net_radiation(i, j, grid, fields, ::AbstractRadiativeFluxes) = fields.surface_net_radiation[i, j]

"""
Base type for turbulent (latent and sensible) heat flux parameterizations.
"""
abstract type AbstractTurbulentFluxes{NF} <: AbstractProcess{NF} end

"""
    sensible_heat_flux(i, j, grid, fields, ::AbstractTurbulentFluxes)

Return the current sensible heat flux at the given indices.
"""
@propagate_inbounds sensible_heat_flux(i, j, grid, fields, ::AbstractTurbulentFluxes) = fields.sensible_heat_flux[i, j]

"""
    latent_heat_flux(i, j, grid, fields, ::AbstractTurbulentFluxes)

Return the current latent heat flux at the given indices.
"""
@propagate_inbounds latent_heat_flux(i, j, grid, fields, ::AbstractTurbulentFluxes) = fields.latent_heat_flux[i, j]

# Surface energy balance parameterizations

"""
Base type for surface albedo and emissivity parameterizations.
"""
abstract type AbstractAlbedo{NF} end

"""
    albedo(i, j, grid, fields, ::AbstractAlbedo)

Return the current albedo at the given indices.
"""
@propagate_inbounds albedo(i, j, grid, fields, ::AbstractAlbedo) = fields.albedo[i, j]

"""
    emissivity(i, j, grid, fields, ::AbstractAlbedo)

Return the current emissivity at the given indices.
"""
@propagate_inbounds emissivity(i, j, grid, fields, ::AbstractAlbedo) = fields.emissivity[i, j]
