"""
    $TYPEDEF

Standard implementation of the surface energy balance (SEB) that computes the radiative,
turbulent, and ground energy fluxes at the surface. The SEB is also responsible for defining
and solving the so-called *skin temperature* (effective emission temperature of the land surface)
as well as the albedo.
"""
struct SurfaceEnergyBalance{
        NF,
        SkinTemperature <: AbstractSkinTemperature{NF},
        TurbulentFluxes <: AbstractTurbulentFluxes{NF},
        RadiativeFluxes <: AbstractRadiativeFluxes{NF},
        Albedo <: AbstractAlbedo{NF},
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

""" $TYPEDSIGNATURES """
@inline function compute_auxiliary!(
        state, grid,
        seb::SurfaceEnergyBalance,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        hydrology::Optional{AbstractSurfaceHydrology} = nothing,
        snow::Optional{AbstractSnow} = nothing,
        args...
    )
    # diagnose the (optionally snow-aware) albedo, then solve the surface energy balance
    compute_auxiliary!(state, grid, seb.albedo, snow)
    solve_surface_energy_balance!(state, grid, seb, constants, atmos, hydrology, snow)
    return nothing
end

""" $TYPEDSIGNATURES """
initialize!(state, grid, seb::SurfaceEnergyBalance, args...) = initialize!(state, grid, seb.skin_temperature, args...)

"""
    $TYPEDSIGNATURES

Solve the surface energy balance for skin temperature on `grid` based on the current atmospheric
and surface hydrology state.
"""
function solve_surface_energy_balance!(
        state, grid,
        seb::SurfaceEnergyBalance{NF, <:ImplicitSkinTemperature},
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        hydrology::Optional{AbstractSurfaceHydrology} = nothing,
        snow::Optional{AbstractSnow} = nothing,
        args...
    ) where {NF}
    evtr = isnothing(hydrology) ? nothing : get_evapotranspiration(hydrology)
    # Construct outputs as auxiliaries + skin temperature (which is prognostic)
    out = (skin_temperature = state.skin_temperature, auxiliary_fields(state, seb)...)
    # Merge the snow thermal auxiliaries so the (optionally snow-aware) conduction target can be evaluated
    fields = merge(get_fields(state, seb, atmos, evtr), get_fields(state, snow))
    launch!(grid, XY, solve_surface_energy_balance_kernel!, out, fields, seb, constants, atmos, evtr, snow, args...)
    return nothing
end

"""
    $TYPEDSIGNATURES

Compute the surface energy fluxes on `grid` based on the current atmospheric state.
"""
function compute_surface_energy_fluxes!(
        state, grid,
        seb::SurfaceEnergyBalance,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        hydrology::Optional{AbstractSurfaceHydrology} = nothing,
        args...
    )
    evtr = isnothing(hydrology) ? nothing : get_evapotranspiration(hydrology)
    # Construct outputs as auxiliaries + skin temperature (which is prognostic)
    out = (skin_temperature = state.skin_temperature, auxiliary_fields(state, seb)...)
    fields = get_fields(state, seb, atmos, evtr)
    launch!(grid, XY, compute_surface_energy_fluxes_kernel!, out, fields, seb, constants, atmos, evtr, args...)
    return nothing
end

# Kernel functions

"""
    $TYPEDSIGNATURES

Fused kernel function that computes the radiative and turbulent fluxes, as well as the ground heat flux based on the current
skin temperature and humidity fluxes.
"""
@propagate_inbounds function compute_surface_energy_fluxes!(
        out, i, j, grid, fields,
        seb::SurfaceEnergyBalance,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        evtr::Optional{AbstractEvapotranspiration} = nothing,
        snow::Optional{AbstractSnow} = nothing,
        args...
    )
    # Compute radiative fluxes
    radiative_fluxes = compute_surface_upwelling_radiation(i, j, grid, fields, seb.radiative_fluxes, seb.skin_temperature, seb.albedo, constants, atmos)
    out.surface_shortwave_up[i, j, 1] = radiative_fluxes.surface_shortwave_up
    out.surface_longwave_up[i, j, 1] = radiative_fluxes.surface_longwave_up
    out.surface_net_radiation[i, j, 1] = compute_surface_net_radiation(i, j, grid, fields, seb.radiative_fluxes, atmos)
    # Compute turbulent fluxes; `snow` partitions the latent flux (evaporation vs. sublimation) by area fraction
    out.sensible_heat_flux[i, j, 1] = compute_sensible_heat_flux(i, j, grid, fields, seb.turbulent_fluxes, seb.skin_temperature, constants, atmos)
    out.latent_heat_flux[i, j, 1] = compute_latent_heat_flux(i, j, grid, fields, seb.turbulent_fluxes, seb.skin_temperature, constants, atmos, evtr, snow)
    # Compute ground heat flux
    out.ground_heat_flux[i, j, 1] = compute_ground_heat_flux(i, j, grid, fields, seb.skin_temperature, seb)
    return nothing
end

# Kernels (fused)

@kernel inbounds = true function solve_surface_energy_balance_kernel!(
        out, grid, fields,
        seb::SurfaceEnergyBalance,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        evtr::Optional{AbstractEvapotranspiration} = nothing,
        snow::Optional{AbstractSnow} = nothing,
        args...
    )
    i, j = @index(Global, NTuple)

    # Solve for skin temperature; `snow` (after `seb`) makes the conduction target and latent flux snow-aware
    solve_skin_temperature!(out, i, j, grid, fields, seb.skin_temperature, seb, snow, constants, atmos, args...)
    if !isnothing(evtr)
        # Recompute evapotranspiration component fluxes from final skin temperature; `snow` scales the
        # ground evaporation by the snow-free fraction (1 − f_snow)
        out_ET = auxiliary_fields(fields, evtr)
        compute_evapotranspiration_fluxes!(out_ET, i, j, grid, fields, evtr, constants, atmos, snow, args...)
    end
    # Recompute fluxes from final skin temperature; `snow` partitions the latent flux (evaporation vs. sublimation)
    compute_surface_energy_fluxes!(out, i, j, grid, fields, seb, constants, atmos, nothing, snow, args...)
end

@kernel inbounds = true function compute_surface_energy_fluxes_kernel!(
        out, grid, fields,
        seb::SurfaceEnergyBalance,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        args...
    )
    i, j = @index(Global, NTuple)

    # Compute fluxes based on current skin temperature
    compute_surface_energy_fluxes!(out, i, j, grid, fields, seb, constants, atmos, args...)
end
