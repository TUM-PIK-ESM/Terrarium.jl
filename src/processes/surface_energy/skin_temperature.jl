"""
    $TYPEDEF

Simple scheme for prescribed skin temperatures from input variables.
"""
struct PrescribedSkinTemperature <: AbstractSkinTemperature end

variables(::PrescribedSkinTemperature) = (
    auxiliary(:ground_heat_flux, XY(), units=u"W/m^2", desc="Ground heat flux"),
    input(:skin_temperature, XY(), units=u"°C", desc="Longwave emission temperature of the land surface in °C"),
)

"""
    $TYPEDEF

Represents prognostic skin temperature evolved according to:

```math
c Δz \\frac{\\partial T_s}{\\partial t} = R_{\\text{net}} + H_s + H_l + G
```
where ``R_{\\text{net}}`` is the net radiation budget, ``H_s`` is the sensible heat flux, ``H_l`` is the latent
heat flux from sublimation and evapotranspiration, and ``G`` is the ground heat flux.
"""
@kwdef struct PrognosticSkinTemperature{NF} <: AbstractSkinTemperature
    "Effective thermal conductivity at the surface [W m⁻¹ K⁻¹]"
    κₛ::NF = 2.0
end

PrognosticSkinTemperature(::Type{NF}; kwargs...) where {NF} = PrognosticSkinTemperature{NF}(; kwargs...)

variables(::PrognosticSkinTemperature) = (
    prognostic(:skin_temperature, XY(), units=u"°C", desc="Longwave emission temperature of the land surface in °C"),
    auxiliary(:ground_heat_flux, XY(), units=u"W/m^2", desc="Ground heat flux"),
    input(:ground_temperature, XY(), units=u"°C", desc="Temperature of the uppermost ground or soil grid cell in °C")
)

function update_skin_temperature!(state, model::AbstractSurfaceEnergyBalanceModel, skinT::PrognosticSkinTemperature)
    grid = get_grid(model)
    launch!(grid, :xy, update_skin_temperature_kernel!, state, grid, skinT)   
end

"""
Default `compute_auxiliary!` for all skin temperature implementation that computes the ground heat flux from the
current radiative, sensible, and latent fluxes.
"""
function compute_auxiliary!(state, model::AbstractSurfaceEnergyBalanceModel, skinT::AbstractSkinTemperature)
    grid = get_grid(model)
    launch!(grid, :xy, compute_ground_heat_flux!, state, grid, skinT)
end

# Kernels

"""
    compute_ground_heat_flux!(state, grid, ::AbstractSkinTemperature)

Diagnose the ground heat flux as the residual of the net radiation budget and turbulent fluxes.
"""
@kernel function compute_ground_heat_flux!(state, grid, ::AbstractSkinTemperature)
    i, j = @index(Global, NTuple)
    # compute flux terms
    R_net = state.RadNet[i, j]
    H_s = state.sensible_heat_flux[i, j]
    H_l = state.latent_heat_flux[i, j]
    # compute ground heat flux as residual of R_net and turbulent fluxes
    state.ground_heat_flux[i, j, 1] = R_net - H_s - H_l
end

"""
    update_skin_temperature_kernel!(state, grid, skinT::PrognosticSkinTemperature)

Diagnose the skin temperature implied by the current `ground_heat_flux` and `ground_temperature`.
"""
@kernel function update_skin_temperature_kernel!(state, grid, skinT::PrognosticSkinTemperature)
    i, j = @index(Global, NTuple)
    idx = (i, j)
    # get thickness of topmost soil/ground grid cell
    field_grid = get_field_grid(grid)
    Δz₁ = Δzᵃᵃᶜ(i, j, field_grid.Nz, field_grid)
    # get ground heat flux
    G = state.ground_heat_flux[i, j]
    # get ground temperature
    Tₛ = state.ground_temperature[i, j]
    # compute new skin temperature T₀ by setting G equal to the half-cell heat flux
    κₛ = skin_thermal_conductivity(idx, state, skinT)
    state.skin_temperature[i, j, 1] = Tₛ - G * Δz₁ / (2*κₛ)
end

# Kernel functions

@inline skin_thermal_conductivity(idx, state, skinT::PrognosticSkinTemperature) = skinT.κₛ

@inline @inbounds skin_temperature(idx, state, ::AbstractSkinTemperature) = state.skin_temperature[idx]

@inline @inbounds ground_heat_flux(idx, state, ::AbstractSkinTemperature) = state.ground_heat_flux[idx]
