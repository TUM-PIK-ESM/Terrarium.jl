abstract type AbstractHeatOperator end
struct VerticalHeatConduction <: AbstractHeatOperator end

@kwdef struct SoilThermalProperties{NF}
    κh_w::NF = 0.57 # thermal conductivity of water [Hillel (1982)]
    κh_i::NF = 2.2 # thermal conductivity of ice [Hillel (1982)]
    κh_a::NF = 0.025 # thermal conductivity of air [Hillel (1982)]
    κh_m::NF = 3.8 # thermal conductivity of mineral soil [Hillel (1982)]
    κh_o::NF = 0.25 # thermal conductivity of organic soil [Hillel (1982)]
    ch_w::NF = 4.2e6 # heat capacity of water
    ch_i::NF = 1.9e6 # heat capacity of ice
    ch_a::NF = 0.00125e6 # heat capacity of air
    ch_m::NF = 2.0e6 # heat capacity of mineral soil
    ch_o::NF = 2.5e6 # heat capacity of organic soil
end

@kwdef struct SoilEnergyBalance{
    HeatOperator<:AbstractHeatOperator
} <: AbstractSoilEnergyBalance
    "Heat transport operator"
    operator::HeatOperator = VerticalHeatConduction()

    "Soil thermal properties"
    thermal_properties::SoilThermalProperties = SoilThermalProperties()
end

@inline function thermalconductivity(i, j, k, state, model::AbstractSoilModel, energy::SoilEnergyBalance)
    # TODO: compute based on soil stratigraphy
    return one(model.grid)
end

@inline function heatcapacity(i, j, k, state, model::AbstractSoilModel, energy::SoilEnergyBalance)
    # TODO: compute based on soil stratigraphy
    return one(model.grid)
end

@inline function update_state(i, j, k, state, model, energy::SoilEnergyBalance)
    # TODO: freeze/thaw
end

@inline function compute_tendencies(i, j, k, state, model, energy::SoilEnergyBalance)
    # Get thermal conductivity
    κ = thermalconductivity(i, j, k, state, model, energy)
    
    # Interior heat flux
    dHdt = -∂zᵃᵃᶜ(i, j, k, grid, diffusive_flux, κ, T)

    # Ground heat flux (upper boundary)
    Qg = @inbounds state.ground_heat_flux[i, j, k]
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    dHdt += ifelse(k == grid.Nz, -Qg / Δz, zero(grid))

    # Geothermal heat flux (lower boundary)
    # TODO: should this just be a flux boundary condition on the tendency?
    Qgeo = @inbounds state.geothermal_heat_flux[i, j, k]
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    dHdt += ifelse(k == 1, Qgeo / Δz, zero(grid))
end

@inline diffusive_flux(i, j, k, grid, κ, T) = -κ * ∂zᵃᵃᶠ(i, j, k, grid, T)
