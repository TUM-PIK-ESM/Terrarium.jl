abstract type AbstractEvapotranspirativeResistance end

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

sensible_heat_flux(idx, state, ::PrescribedTurbulentFluxes) = state.sensible_heat_flux[idx]

latent_heat_flux(idx, state, ::PrescribedTurbulentFluxes) = state.latent_heat_flux[idx]

"""
    $TYPEDEF

Represents the standard case where the turbulent (sensible and latent) heat fluxes are diagnosed
from atmosphere and soil conditions.
"""
struct DiagnosedTurbulentFluxes{
    ER<:AbstractEvapotranspirativeResistance,
    LF<:AbstractSoilMoistureLimitingFactor
} <: AbstractTurbulentFluxes
    "Parameterization of aerodynamic (rₐ) and canopy/surface (rₛ) resistances"
    resistance::ER

    "Parameterization of the soil moisture limiting factor [-]"
    soil_moisture_limiting_factor::LF
end

function DiagnosedTurbulentFluxes(
    ::Type{NF};
    resistance::AbstractEvapotranspirativeResistance=ConstantAerodynamicResistance(NF),
    soil_moisture_limiting_factor::AbstractSoilMoistureLimitingFactor=UnlimitedSoilMoisture()
) where {NF}
    return DiagnosedTurbulentFluxes(resistance, soil_moisture_limiting_factor)
end

variables(tur::DiagnosedTurbulentFluxes) = (
    auxiliary(:sensible_heat_flux, XY(), units=u"W/m^2", desc="Sensible heat flux at the surface [W m⁻²]"),
    auxiliary(:latent_heat_flux, XY(), units=u"W/m^2", desc="Latent heat flux at the surface [W m⁻²]"),
    variables(tur.resistance)...,
    variables(tur.soil_moisture_limiting_factor)...,
)

function compute_auxiliary!(state, model, tur::DiagnosedTurbulentFluxes)
    (; grid, surface_energy_balance, atmosphere, constants) = model
    launch!(
        state,
        grid,
        :xy,
        compute_turbulent_fluxes!,
        tur,
        surface_energy_balance.skin_temperature,
        atmosphere,
        constants,
    )
end

@kernel function compute_turbulent_fluxes!(
    state,
    ::AbstractLandGrid,
    tur::DiagnosedTurbulentFluxes,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    i, j = @index(Global, NTuple)
    idx = (i, j)
    # compute sensible heat flux
    state.sensible_heat_flux[i, j, 1] = sensible_heat_flux(idx, state, tur, skinT, atmos, constants)
    # compute latent heat flux
    state.latent_heat_flux[i, j, 1] = latent_heat_flux(idx, state, tur, skinT, atmos, constants)
end

# Kernel functions

function sensible_heat_flux(
    idx, state,
    tur::DiagnosedTurbulentFluxes,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    let ρₐ = constants.ρₐ, # density of air
        cₐ = constants.cₐ, # specific heat capacity of air
        rₐ = aerodynamic_resistance(idx, state, tur.resistance), # aerodynamic resistance
        Ts = skin_temperature(idx, state, skinT), # skin temperature
        Tair = air_temperature(idx, state, atmos); # air temperature
        # Calculate sensible heat flux (positive upwards)
        Hₛ = -ρₐ * cₐ / rₐ * (Tair - Ts)
        return Hₛ
    end
end

function latent_heat_flux(
    idx, state,
    tur::DiagnosedTurbulentFluxes,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    let L = constants.Lsg, # specific latent heat of vaporization of water
        ρₐ = constants.ρₐ, # density of air
        rₐ = aerodynamic_resistance(idx, state, tur.resistance), # aerodynamic resistance
        β = soil_moisture_limiting_factor(idx, state, tur.soil_moisture_limiting_factor),
        q_sat = surface_humidity_at_saturation(idx, state, skinT, atmos, constants), # near-surface specific humidity
        q_air = specific_humidity(idx, state, atmos); # specific humidity of the air
        # Calculate latent heat flux (positive upwards)
        Hₗ = -L * ρₐ * β / rₐ * (q_air - q_sat)
        return Hₗ
    end
end

"""
Near-surface specific humidity at saturation, ``q_{\\text{sat}}``, determined based on
the skin temperature and atmospheric pressure.
"""
function surface_humidity_at_saturation(
    idx, state,
    skinT::AbstractSkinTemperature,
    atmos::AbstractAtmosphere,
    constants::PhysicalConstants
)
    let T = skin_temperature(idx, state, skinT),
        γ = constants.γ,
        p = air_pressure(idx, state, atmos),
        e = saturation_vapor_pressure(T);
        # convert saturation vapor pressure to specific humidity
        q_sat = γ * e / p
        return q_sat
    end
end

"""
    $TYPEDEF

Dummy implementation of the aerodynamic resistance that simply returns a constant value.
"""
@kwdef struct ConstantAerodynamicResistance{NF} <: AbstractEvapotranspirativeResistance
    "Constant aerodynamic resistance [s m⁻¹]"
    rₐ::NF = 50.0
end

ConstantAerodynamicResistance(::Type{NF}; kwargs...) where {NF} = ConstantAerodynamicResistance{NF}(; kwargs...)

aerodynamic_resistance(idx, state, res::ConstantAerodynamicResistance) = res.rₐ

struct UnlimitedSoilMoisture <: AbstractSoilMoistureLimitingFactor end

soil_moisture_limiting_factor(idx, state, ::UnlimitedSoilMoisture) = 1
