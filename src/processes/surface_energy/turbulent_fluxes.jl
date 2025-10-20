abstract type AbstractEvapotranspirativeResistance end

abstract type AbstractSurfaceHumidity end

abstract type AbstractSoilMoistureLimitingFactor end

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

sensible_heat_flux(idx, state, ::PrescribedTurbulentFluxes) = state.sensible_heat_flux[idx...]

latent_heat_flux(idx, state, ::PrescribedTurbulentFluxes) = state.latent_heat_flux[idx...]

"""
    $TYPEDEF

Represents the standard case where the turbulent (sensible and latent) heat fluxes are diagnosed
from atmosphere and soil conditions.
"""
struct DiagnosedTurbulentFluxes{
    ER<:AbstractEvapotranspirativeResistance,
    SH<:AbstractSurfaceHumidity,
    LF<:AbstractSoilMoistureLimitingFactor
} <: AbstractTurbulentFluxes
    "Parameterization of aerodynamic (rₐ) and canopy/surface (rₛ) resistances"
    resistance::ER

    "Parameterization of surface humidity (qₛ) [kg kg⁻¹]"
    surface_humidity::SH

    "Parameterization of the soil moisture limiting factor [-]"
    soil_moisture_limiting_factor::LF
end

function DiagnosedTurbulentFluxes(
    ::Type{NF};
    resistance::AbstractEvapotranspirativeResistance=AerodynamicResistance(NF),
    surface_humidity::AbstractSurfaceHumidity=SurfaceHumidityWaterIce(NF),
    soil_moisture_limiting_factor::AbstractSoilMoistureLimitingFactor=UnlimitedSoilMoisture()
) where {NF}
    return DiagnosedTurbulentFluxes(resistance, surface_humidity, soil_moisture_limiting_factor)
end

variables(tur::DiagnosedTurbulentFluxes) = (
    auxiliary(:sensible_heat_flux, XY(), units=u"W/m^2", desc="Sensible heat flux at the surface [W m⁻²]"),
    auxiliary(:latent_heat_flux, XY(), units=u"W/m^2", desc="Latent heat flux at the surface [W m⁻²]"),
    variables(tur.resistance)...,
    variables(tur.surface_humidity)...,
    variables(tur.soil_moisture_limiting_factor)...,
)

function compute_auxiiliary!(state, model, tur::DiagnosedTurbulentFluxes)
    grid = get_grid(model)
    launch!(grid, compute_turbulent_fluxes!, state, grid, tur)
end

@kernel function compute_turbulent_fluxes!(state, grid, fluxes::DiagnosedTurbulentFluxes)
    i, j = idx = @index(Global, NTuple)
    # compute sensible heat flux
    state.sensible_heat_flux[i, j, 1] = sensible_heat_flux(idx, state, fluxes)
    # compute latent heat flux
    state.latent_heat_flux[i, j, 1] = latent_heat_flux(idx, state, fluxes)
end

# Kernel functions

function sensible_heat_flux(
    idx, state,
    fluxes::DiagnosedTurbulentFluxes,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    
end

function latent_heat_flux(
    idx, state,
    fluxes::DiagnosedTurbulentFluxes,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)

end

@kwdef struct ConstantAerodynamicResistance{NF} <: AbstractEvapotranspirativeResistance
    "Constant aerodynamic resistance [s m⁻¹]"
    rₐ::NF = 50.0
end

aerodynamic_resistance(idx, state, res::ConstantAerodynamicResistance) = res.rₐ

struct SurfaceHumidityWaterIce <: AbstractSurfaceHumidity end

"""
Surface humidity at saturation ``q_{\\text{sat}}``.
"""
function surface_humidity_at_saturation(
    idx, state,
    ::SurfaceHumidityWaterIce,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    let T = skin_temperature(idx, state, skinT),
        γ = constants.γ,
        p = air_pressure(idx, state, atmos);
        q_sat = γ * saturation_vapor_pressure(T) / p
        return q_sat
    end
end
