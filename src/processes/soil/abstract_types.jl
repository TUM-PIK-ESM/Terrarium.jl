# Component types

"""
    $TYPEDEF

Base type for coupled ground processes.
"""
abstract type AbstractGround{NF} <: AbstractCoupledProcesses{NF} end

"""
    get_stratigraphy(ground)

Return the ground stratigraphy parameterization associated with `ground`.
"""
@inline get_stratigraphy(ground::AbstractGround) = ground.strat

"""
    get_energy(ground)

Return the energy balance scheme associated with `ground`.
"""
@inline get_energy(ground::AbstractGround) = ground.energy

"""
    get_hydrology(ground)

Return the hydrology scheme associated with `ground`.
"""
@inline get_hydrology(ground::AbstractGround) = ground.hydrology

"""
    $TYPEDEF

Base type for coupled soil processes.
"""
abstract type AbstractSoil{NF} <: AbstractGround{NF} end

"""
    get_biogeochemistry(soil)

Return the biogeochemistry scheme associated with `soil`.
"""
@inline get_biogeochemistry(soil::AbstractSoil) = soil.biogeochem

# Soil process types

"""
    $TYPEDEF

Base type for soil energy balance process implementations. Subtypes should define
state variables for soil `temperature`, `internal_energy`, and any other relevant
thermal properties or state variables.
"""
abstract type AbstractSoilEnergyBalance{NF} <: AbstractProcess{NF} end

"""
    compute_energy_tendency(i, j, k, grid, ::SoilEnergyBalance, args...)

Compute the internal energy tendency `∂U∂t` at index `i, j, k`.
"""
function compute_energy_tendency end

"""
    compute_thermal_conductivity(i, j, k, grid, ::SoilEnergyBalance, args...)

Compute the thermal conductivity at index `i, j, k`.
"""
function compute_thermal_conductivity end

"""
    $TYPEDEF

Base type for soil hydrology implementations. Subtypes should define
state variables for `saturation_water_ice`, `hydraulic_conductivity`, 
`liquid_water_fraction`, and the current `water_table` level, along with
any other implementation-specific state variables.
"""
abstract type AbstractSoilHydrology{NF} <: AbstractProcess{NF} end

"""
    get_swrc(hydrology::AbstractSoilHydrology)

Return the soil water retention curve from the `hydraulic_properties` associated with
the given soil hydrology configuration.
"""
function get_swrc end

"""
    get_hydraulic_properties(hydrology::AbstractSoilHydrology)

Return the soil hydraulic properties defined by the given soil `hydrology` configuration.
"""
function get_hydraulic_properties end

"""
    saturation_water_ice(i, j, k, grid, fields, ::AbstractSoilHydrology)

Compute or retrieve the current saturation level of water + ice in the pore space.
"""
function saturation_water_ice end

"""
    hydraulic_conductivity(i, j, k, grid, fields, ::AbstractSoilHydrology)

Compute or retrieve the current hydraulic conductivity at grid cell `i, j` and vertical
layer **face** `k`.
"""
function hydraulic_conductivity end

"""
    liquid_water_fraction(i, j, k, grid, fields, ::AbstractSoilHydrology)

Compute or retrieve the current fraction of unfrozen water in the pore space.
"""
function liquid_water_fraction end

"""
    water_table(i, j, k, grid, fields, ::AbstractSoilHydrology)

Compute or retrieve the current water table level relative to the surface.
"""
function water_table end

"""
    surface_excess_water(i, j, k, grid, fields, ::AbstractSoilHydrology)

Retrieve the current saturation level of water + ice in the pore space.
"""
function surface_excess_water end

"""
    $TYPEDEF

Base type for soil biogeochemistry process implementations.
"""
abstract type AbstractSoilBiogeochemistry{NF} <: AbstractProcess{NF} end

"""
    density_soc(i, j, k, grid, fields, bgc::AbstractSoilBiogeochemistry)

Compute or return the soil organic carbon density at grid cell index `i, j, k`.
"""
function density_soc end

# Parameterization types

"""
    $TYPEDEF

Base type for mineral soil texture parameterizations.
"""
abstract type AbstractSoilTexture{NF} end

"""
    $TYPEDEF

Base type for parameterizations of soil porosity.
"""
abstract type AbstractSoilPorosity{NF} end

"""
    mineral_porosity(::AbstractSoilPorosity, texture::SoilTexture)

Compute or retrieve the natural porosity of the mineral soil constitutents, i.e.
excluding organic material.
"""
function mineral_porosity end

"""
    organic_porosity(::AbstractSoilPorosity, texture::SoilTexture)

Compute or retrieve the natural porosity of the organic soil constitutents, i.e.
excluding mineral material.
"""
function organic_porosity end

"""
    $TYPEDEF

Base type for soil stratigraphy parameterizations.
"""
abstract type AbstractStratigraphy{NF} end

"""
    soil_texture(i, j, k, grid, fields, ::AbstractStratigraphy, args...)

Return the texture of the soil at index `i, j, k` for the given stratigraphy parameterization.
"""
function soil_texture end

"""
    soil_matrix(i, j, k, grid, fields, ::AbstractStratigraphy, args...)

Return the solid matrix of the soil at index `i, j, k` for the given stratigraphy parameterization.
"""
function soil_matrix end

"""
    soil_volume(i, j, k, grid, fields, ::AbstractStratigraphy, args...)

Return a [`SoilVolume`](@ref) describing the full material composition of the soil volume at index
`i, j, k` for the given stratigraphy parameterization.
"""
function soil_volume end

"""
    $TYPEDEF

Base type for formulations of the heat transfer operator.
"""
abstract type AbstractHeatOperator end

"""
    $TYPEDEF

Base type for closure relations between energy and temperature in soil volumes.
"""
abstract type AbstractSoilEnergyClosure <: AbstractClosureRelation end

"""
    $TYPEDEF

Base type for closure relations between water saturation and potential in soil volumes.
"""
abstract type AbstractSoilWaterClosure <: AbstractClosureRelation end
