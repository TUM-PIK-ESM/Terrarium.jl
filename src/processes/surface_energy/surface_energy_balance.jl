struct SurfaceEnergyBalance{
    Albedo<:AbstractAlbedo,
    SkinTemperature<:AbstractSkinTemperature,
    TurbulentFluxes<:AbstractTurbulentFluxes,
    RadiativeFluxes<:AbstractRadiativeFluxes
} <: AbstractSurfaceEnergyBalance
    "Scheme for parameterizing surface albedo"
    albedo::Albedo

    "Scheme for determining skin temperature and ground heat flux"
    skin_temperature::SkinTemperature

    "Scheme for determining the net radiation budget"
    radiative_fluxes::RadiativeFluxes

    "Scheme for computing turbulent (sensible and latent) heat fluxes"
    turbulent_fluxes::TurbulentFluxes
end

function SurfaceEnergyBalance(
    ::Type{NF};
    albedo::AbstractAlbedo = ConstantAlbedo(NF),
    radiative_fluxes::AbstractRadiativeFluxes = DiagnosedRadiativeFluxes(),
    turbulent_fluxes::AbstractTurbulentFluxes = DiagnosedTurbulentFluxes(NF),
    skin_temperature::AbstractSkinTemperature = ImplicitSkinTemperature()
) where {NF}
    return SurfaceEnergyBalance(albedo, skin_temperature, radiative_fluxes, turbulent_fluxes)
end

processes(seb::SurfaceEnergyBalance) = (
    seb.albedo,
    seb.skin_temperature,
    seb.radiative_fluxes,
    seb.turbulent_fluxes
)

function compute_auxiliary!(state, model, seb::SurfaceEnergyBalance)
    compute_auxiliary!(state, model, seb.albedo)
    compute_auxiliary!(state, model, seb.skin_temperature)
    # compute energy fluxes
    compute_surface_energy_fluxes!(state, model, seb)
end

function compute_auxiliary!(state, model, seb::SurfaceEnergyBalance{AL, <:ImplicitSkinTemperature}) where {AL<:AbstractAlbedo}
    compute_auxiliary!(state, model, seb.albedo)
    # set skin temeprature initially to ground temperature
    set!(state.skin_temperature, state.ground_temperature)
    # compute fluxes
    compute_surface_energy_fluxes!(state, model, seb)
    # diagnose skin temperature
    update_skin_temperature!(state, model, seb.skin_temperature)
    # recompute fluxes
    compute_surface_energy_fluxes!(state, model, seb)
end

function compute_surface_energy_fluxes!(state, model, seb::SurfaceEnergyBalance)
    compute_auxiliary!(state, model, seb.radiative_fluxes)
    compute_auxiliary!(state, model, seb.turbulent_fluxes)
    compute_auxiliary!(state, model, seb.skin_temperature)
end
