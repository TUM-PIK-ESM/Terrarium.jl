# Interface for differential operators
abstract type AbstractOperator end

"""
    get_closure(op::AbstractOperator)

Returns an `AbstractClosureRelation` for the given differential operator.
Deefaults to returning `nothing` (i.e. no closure).
"""
get_closure(op::AbstractOperator)::AbstractClosureRelation = nothing

variables(op::AbstractOperator) = ()

# Interface for processes

"""
    AbstractProcess

Base type for all processes which define physics or parameterizations but are not standalone models.
"""
abstract type AbstractProcess end

variables(process::AbstractProcess) = ()

initialize!(state, model, process::AbstractProcess) = compute_auxiliary!(state, model, process)

compute_auxiliary!(state, model, ::AbstractProcess) = nothing

compute_tendencies!(state, model, ::AbstractProcess) = nothing

# Soil process types
# TODO: Think more about these process types and what methods they should have.

abstract type AbstractStratigraphy{NF} <: AbstractProcess end

abstract type AbstractSoilEnergyBalance{NF} <: AbstractProcess end

abstract type AbstractSoilHydrology{NF} <: AbstractProcess end

abstract type AbstractSoilBiogeochemistry{NF} <: AbstractProcess end

# Vegetation process types

abstract type AbstractPhotosynthesis end

abstract type AbstractStomatalConductance end

abstract type AbstractAutotrophicRespiration end

abstract type AbstractVegetationDynamics end

abstract type AbstractPhenology end

abstract type AbstractVegetationCarbonDynamics end

# Surface energy balance types

"""
Base type for surface albedo and emissivity process implementations.
"""
abstract type AbstractAlbedo <: AbstractProcess end

"""
    alebdo(idx, state, ::AbstractAlbedo)

Return the current albedo at the given `idx`.
"""
function albedo end

"""
    emissivity(idx, state, ::AbstractAlbedo)
AbstractAlbedo
Return the current emissivity at the given `idx`.
"""
function emissivity end

"""
Base type for radiation budget schemes.
"""
abstract type AbstractRadiativeFluxes <: AbstractProcess end

"""
    net_radiation(idx, state, ::AbstractRadiativeFluxes)

Return the current net radiation at the given `idx`.
"""
function net_radiation end

"""
Base type for turbulent (latent and sensible) heat fluxes at the surface.
"""
abstract type AbstractTurbulentFluxes <: AbstractProcess end

"""
    sensible_heat_flux(idx, state, ::AbstractTurbulentFluxes)

Return the current sensible heat flux at the given `idx`.
"""
function sensible_heat_flux end

"""
    latent_heat_flux(idx, state, ::AbstractTurbulentFluxes)

Return the current latent heat flux at the given `idx`.
"""
function latent_heat_flux end

"""
Base type for skin temperature and ground heat flux schemes.
"""
abstract type AbstractSkinTemperature <: AbstractProcess end

"""
    skin_temperature(idx, state, ::AbstractSkinTemperature)

Return the current skin temperature at the given `idx`.
"""
function skin_temperature end

"""
Base type for surface energy balance schemes.
"""
abstract type AbstractSurfaceEnergyBalance <: AbstractProcess end

# Atmosphere

abstract type AbstractHumidity end

abstract type AbstractPrecipitation end

abstract type AbstractIncomingRadiation end

abstract type AbstractAtmosphere{
    PR<:AbstractPrecipitation,
    IR<:AbstractIncomingRadiation,
    HM<:AbstractHumidity
} <: AbstractBoundaryConditions
end
