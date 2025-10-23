"""
    $TYPEDEF

Represents skin temperature prescribed by an input variable.
"""
struct PrescribedSkinTemperature <: AbstractSkinTemperature end

variables(::PrescribedSkinTemperature) = (
    auxiliary(:ground_heat_flux, XY(), units=u"W/m^2", desc="Ground heat flux"),
    input(:skin_temperature, XY(), units=u"°C", desc="Longwave emission temperature of the land surface in °C"),
)

function compute_auxiiliary!(state, model, skinT::AbstractSkinTemperature)
    grid = get_grid(model)
    rad = get_radiative_fluxes(model)
    tur = get_turbulent_fluxes(model)
    launch!(grid, compute_ground_heat_flux!, state, grid, skinT, rad, tur)
end

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

# Kernels

@kernel function compute_ground_heat_flux!(
    state, grid,
    ::PrognosticSkinTemperature,
    rad::AbstractRadiativeFluxes,
    tur::AbstractTurbulentFluxes,
)
    # compute flux terms
    R_net = net_radiation(idx, state, rad)
    H_s = sensible_heat_flux(idx, state, tur)
    H_l = latent_heat_flux(idx, state, tur)
    # compute ground heat flux as residual of R_net and turbulent fluxes
    state.ground_heat_flux[i, j] = -R_net - H_s - H_l
end

"""
    update_skin_temperature!(idx, state, grid, skinT::PrognosticSkinTemperature)

Diagnose the skin temperature implied by the current `ground_heat_flux` and `ground_temperature`.
"""
@kernel function update_skin_temperature!(state, grid, skinT::PrognosticSkinTemperature)
    i, j = @index(Global, NTuple)
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
