"""
    $TYPEDEF

Standard implementation of the surface energy balance (SEB) that computes the radiative,
turbulent, and ground energy fluxes at the surface. The SEB is also responsible for defining
and solving the so-called *skin temperature* (effective emission temperature of the land surface)
as well as the albedo.
"""
struct SurfaceEnergyBalance{
    NF,
    SkinTemperature<:AbstractSkinTemperature{NF},
    TurbulentFluxes<:AbstractTurbulentFluxes{NF},
    RadiativeFluxes<:AbstractRadiativeFluxes{NF},
    Albedo<:AbstractAlbedo{NF}
} <: AbstractSurfaceEnergyBalance{NF}
    "Scheme for determining skin temperature and ground heat flux"
    skin_temperature::SkinTemperature

    "Scheme for determining the net radiation budget"
    radiative_fluxes::RadiativeFluxes

    "Scheme for computing turbulent (sensible and latent) heat fluxes"
    turbulent_fluxes::TurbulentFluxes

    "Scheme for parameterizing surface albedo"
    albedo::Albedo
end

function SurfaceEnergyBalance(
    ::Type{NF};
    radiative_fluxes::AbstractRadiativeFluxes = DiagnosedRadiativeFluxes(NF),
    turbulent_fluxes::AbstractTurbulentFluxes = DiagnosedTurbulentFluxes(NF),
    skin_temperature::AbstractSkinTemperature = ImplicitSkinTemperature(NF),
    albedo::AbstractAlbedo = ConstantAlbedo(NF)
) where {NF}
    return SurfaceEnergyBalance(skin_temperature, radiative_fluxes, turbulent_fluxes, albedo)
end

variables(seb::SurfaceEnergyBalance) = tuplejoin(
    variables(seb.albedo),
    variables(seb.skin_temperature),
    variables(seb.radiative_fluxes),
    variables(seb.turbulent_fluxes)
)

function compute_auxiliary!(
    state, grid,
    seb::SurfaceEnergyBalance,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    compute_surface_energy_fluxes!(state, grid, seb, atmos, constants)
end

"""
    $TYPEDSIGNATURES

Compute the surface energy fluxes on `grid` based on the current atmospheric state.
"""
function compute_surface_energy_fluxes!(
    state, grid,
    seb::SurfaceEnergyBalance,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants,
    args...
)
    launch!(grid, XY, compute_surface_energy_fluxes_kernel!, state, seb, atmos, constants, args...)
end

# Kernels (fused)

@kernel inbounds=true function compute_surface_energy_fluxes_kernel!(
    state, grid,
    seb::SurfaceEnergyBalance,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants,
    args...
)
    i, j = @index(Global, NTuple)

    # Compute fluxes based on current skin temperature
    compute_surface_energy_fluxes!(state, i, j, grid, seb, atmos, constants, args...)
end

@kernel inbounds=true function compute_surface_energy_fluxes_kernel!(
    state, grid,
    seb::SurfaceEnergyBalance{NF, <:ImplicitSkinTemperature},
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants,
    args...
) where {NF}
    i, j = @index(Global, NTuple)

    # Compute fluxes based on current skin temperature
    compute_surface_energy_fluxes!(state, i, j, grid, seb, atmos, constants, args...)
    # Update skin temperature
    state.skin_temperature[i, j, 1] = compute_skin_temperature(i, j, grid, state, seb.skin_temperature)
    # Recompute fluxes from updated skin temperature
    compute_surface_energy_fluxes!(state, i, j, grid, seb, atmos, constants, args...)
end

# Kernel functions

"""
    $TYPEDSIGNATURES

Mutating kernel function that computes the radiative and turbulent fluxes, as well as the ground heat flux, in
a single kernel operation for efficiency.
"""
@propagate_inbounds function compute_surface_energy_fluxes!(
    state, i, j, grid,
    seb::SurfaceEnergyBalance,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    # Compute radiative fluxes
    radiative_fluxes = compute_surface_upwelling_radiation(i, j, grid, state, seb.radiative_fluxes, seb.skin_temperature, seb.albedo, atmos, constants)
    state.surface_shortwave_up[i, j, 1] = radiative_fluxes.surface_shortwave_up[i, j]
    state.surface_longwave_up[i, j, 1] = radiative_fluxes.surface_longwave_up[i, j]
    state.surface_net_radiation[i, j, 1] = compute_surface_net_radiation(i, j, grid, state, seb.radiative_fluxes)
    # Compute turbulent fluxes
    state.sensible_heat_flux[i, j, 1] = compute_sensible_heat_flux(i, j, grid, state, seb.turbulent_fluxes, seb.skin_temperature, atmos, constants)
    state.latent_heat_flux[i, j, 1] = compute_latent_heat_flux(i, j, grid, state, seb.turbulent_fluxes, seb.skin_temperature, atmos, constants)
    # Compute ground heat flux
    state.ground_heat_flux[i, j, 1] = compute_ground_heat_flux(i, j, grid, state, seb.skin_temperature, seb.radiative_fluxes, seb.turbulent_fluxes)
end
