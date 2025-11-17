# Interface for differential operators
abstract type AbstractOperator end

"""
    get_closure(op::AbstractOperator)

Returns an `AbstractClosureRelation` for the given differential operator.
Deefaults to returning `nothing` (i.e. no closure).
"""
get_closure(op::AbstractOperator) = nothing

variables(op::AbstractOperator) = ()

closure!(state, model::AbstractModel, ::Nothing) = nothing

invclosure!(state, model::AbstractModel, ::Nothing) = nothing

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

# also allow dispatch on nothing
compute_auxiliary!(state, model, ::Nothing) = nothing
compute_tendencies!(state, model, ::Nothing) = nothing

# Soil process types
# TODO: Think more about these process types and what methods they should have.

abstract type AbstractStratigraphy{NF} <: AbstractProcess end

abstract type AbstractSoilEnergyBalance{NF} <: AbstractProcess end

abstract type AbstractSoilHydrology{NF} <: AbstractProcess end

abstract type AbstractSoilBiogeochemistry{NF} <: AbstractProcess end

# Vegetation process types

abstract type AbstractPhotosynthesis <: AbstractProcess  end

abstract type AbstractStomatalConductance <: AbstractProcess  end

abstract type AbstractAutotrophicRespiration <: AbstractProcess  end

abstract type AbstractVegetationDynamics <: AbstractProcess  end

abstract type AbstractPhenology <: AbstractProcess  end

abstract type AbstractVegetationCarbonDynamics <: AbstractProcess  end

# Surface energy balance types

"""
Base type for surface albedo and emissivity process implementations.
"""
abstract type AbstractAlbedo <: AbstractProcess end

"""
    albedo(i, j, state, ::AbstractAlbedo)

Return the current albedo at the given indices.
"""
function albedo end

"""
    emissivity(i, j, state, ::AbstractAlbedo)

Return the current emissivity at the given indices.
"""
function emissivity end

"""
Base type for radiation budget schemes.
"""
abstract type AbstractRadiativeFluxes <: AbstractProcess end

"""
    surface_net_radiation(i, j, state, ::AbstractRadiativeFluxes)

Return the current net radiation at the given indices.
"""
function surface_net_radiation end

"""
Base type for turbulent (latent and sensible) heat fluxes at the surface.
"""
abstract type AbstractTurbulentFluxes <: AbstractProcess end

"""
    sensible_heat_flux(i, j, state, ::AbstractTurbulentFluxes)

Return the current sensible heat flux at the given indices.
"""
function sensible_heat_flux end

"""
    latent_heat_flux(i, j, state, ::AbstractTurbulentFluxes)

Return the current latent heat flux at the given indices.
"""
function latent_heat_flux end

"""
Base type for skin temperature and ground heat flux schemes.
"""
abstract type AbstractSkinTemperature <: AbstractProcess end

"""
    skin_temperature(i, j, state, ::AbstractSkinTemperature)

Return the current skin temperature at the given indices.
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
} <: AbstractProcess
end
