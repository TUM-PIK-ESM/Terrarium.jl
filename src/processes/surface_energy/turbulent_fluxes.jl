"""
    $TYPEDEF

Represents the simplest case where the turbulent (sensible and latent) heat fluxes are prescribed
via input variables.
"""
struct PrescribedTurbulentFluxes <: AbstractTurbulentFluxes end

variables(::PrescribedTurbulentFluxes) = (
    input(:sensible_heat_flux, XY(), units=u"W/m^2", desc="Sensible heat flux at the surface [W m⁻²]"),
    input(:latent_heat_flux, XY(), units=u"W/m^2", desc="Latent heat flux at the surface [W m⁻²]")
)

sensible_heat_flux(i, j, grid, state, ::PrescribedTurbulentFluxes) = state.sensible_heat_flux[i, j]

latent_heat_flux(i, j, grid, state, ::PrescribedTurbulentFluxes) = state.latent_heat_flux[i, j]

"""
    $TYPEDEF

Represents the standard case where the turbulent (sensible and latent) heat fluxes are diagnosed
from atmosphere and soil conditions.
"""
struct DiagnosedTurbulentFluxes{NF} <: AbstractTurbulentFluxes end

DiagnosedTurbulentFluxes(::Type{NF}) where {NF} = DiagnosedTurbulentFluxes{NF}()

# Process methods

variables(::DiagnosedTurbulentFluxes) = (
    auxiliary(:sensible_heat_flux, XY(), units=u"W/m^2", desc="Sensible heat flux at the surface [W m⁻²]"),
    auxiliary(:latent_heat_flux, XY(), units=u"W/m^2", desc="Latent heat flux at the surface [W m⁻²]"),
)

function compute_auxiliary!(state, model, tur::DiagnosedTurbulentFluxes)
    (; grid, surface_energy_balance, atmosphere, constants) = model
    launch!(grid, XY, compute_turbulent_fluxes_kernel!, state,
            tur, surface_energy_balance.skin_temperature, atmosphere, constants)
end

# Kernels

@kernel function compute_turbulent_fluxes_kernel!(
    state, grid,
    tur::DiagnosedTurbulentFluxes,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    i, j = @index(Global, NTuple)
    # compute sensible heat flux
    state.sensible_heat_flux[i, j, 1] = sensible_heat_flux(i, j, grid, state, tur, skinT, atmos, constants)
    # compute latent heat flux
    state.latent_heat_flux[i, j, 1] = latent_heat_flux(i, j, grid, state, tur, skinT, atmos, constants)
end

@kernel function compute_turbulent_fluxes_kernel!(
    state, grid,
    tur::DiagnosedTurbulentFluxes,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants,
    evtr::AbstractEvapotranspiration
)
    i, j = @index(Global, NTuple)
    # compute sensible heat flux
    state.sensible_heat_flux[i, j, 1] = sensible_heat_flux(i, j, grid, state, tur, skinT, atmos, constants)
    # compute latent heat flux
    state.latent_heat_flux[i, j, 1] = latent_heat_flux(i, j, grid, state, tur, evtr)
end

# Kernel functions

@inline function sensible_heat_flux(
    i, j, grid, state,
    tur::DiagnosedTurbulentFluxes,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    let ρₐ = constants.ρₐ, # density of air
        cₐ = constants.cₐ, # specific heat capacity of air
        rₐ = aerodynamic_resistance(i, j, grid, state, atmos), # aerodynamic resistance
        Ts = skin_temperature(i, j, grid, state, skinT), # skin temperature
        Tair = air_temperature(i, j, grid, state, atmos); # air temperature
        # Calculate sensible heat flux (positive upwards)
        Hₛ = -ρₐ * cₐ / rₐ * (Tair - Ts)
        return Hₛ
    end
end

# Standalone latent heat flux
@inline function latent_heat_flux(
    i, j, grid, state,
    ::DiagnosedTurbulentFluxes,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    let L = constants.Llg, # specific latent heat of vaporization of water
        ρₐ = constants.ρₐ, # density of air
        Ts = skin_temperature(i, j, grid, state, skinT),
        rₐ = aerodynamic_resistance(i, j, grid, state, atmos), # aerodynamic resistance
        Δq = compute_humidity_vpd(i, j, grid, state, atmos, constants, Ts);
        # Calculate latent heat flux (positive upwards)
        Hₗ = L * ρₐ * Δq / rₐ
        return Hₗ
    end
end

# Latent heat flux based on ET
@inline function latent_heat_flux(
    i, j, grid, state,
    ::DiagnosedTurbulentFluxes,
    evtr::AbstractEvapotranspiration
)
    let L = constants.Llg, # specific latent heat of vaporization of water
        ρₐ = constants.ρₐ, # density of air
        Qh = surface_humidity_flux(i, j, grid, state, evtr);
        # Calculate latent heat flux (positive upwards)
        Hₗ = L * ρₐ * Qh
        return Hₗ
    end
end
