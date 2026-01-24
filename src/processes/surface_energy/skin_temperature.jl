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

Scheme for an implicit skin temperature ``T_0`` satisfying:
```math
R_{\\text{net}}(T_0) = H_s(T_0) + H_l(T_0) + G(T_0, T_1)
```
where ``R_{\\text{net}}`` is the net radiation budget, ``H_s`` is the sensible heat flux, ``H_l`` is the latent
heat flux from sublimation and evapotranspiration, ``G`` is the ground heat flux, and ``T_1`` is the ground
temperature, or temperature of the uppermost subsurface (soil or snow) layer.
"""
@kwdef struct ImplicitSkinTemperature{NF} <: AbstractSkinTemperature
    "Effective thermal conductivity at the surface [W m⁻¹ K⁻¹]"
    κₛ::NF = 2.0
end

ImplicitSkinTemperature(::Type{NF}; kwargs...) where {NF} = ImplicitSkinTemperature{NF}(; kwargs...)

variables(::ImplicitSkinTemperature) = (
    prognostic(:skin_temperature, XY(), units=u"°C", desc="Longwave emission temperature of the land surface in °C"),
    auxiliary(:ground_heat_flux, XY(), units=u"W/m^2", desc="Ground heat flux"),
    input(:ground_temperature, XY(), units=u"°C", desc="Temperature of the uppermost ground or soil grid cell in °C")
)

"""
Default `compute_auxiliary!` for all skin temperature implementation that computes the ground heat flux from the
current radiative, sensible, and latent fluxes.
"""
function compute_auxiliary!(state, model, skinT::AbstractSkinTemperature)
    compute_ground_heat_flux!(state, model, skinT)
end

function compute_ground_heat_flux!(state, model, skinT::AbstractSkinTemperature)
    grid = get_grid(model)
    launch!(grid, XY, compute_ground_heat_flux_kernel!, state, skinT)
end

function update_skin_temperature!(state, model, skinT::AbstractSkinTemperature)
    grid = get_grid(model)
    launch!(grid, XY, update_skin_temperature_kernel!, state, skinT)
end

# Kernels

"""
    compute_ground_heat_flux_kernel!(state, grid, ::AbstractSkinTemperature)

Diagnose the ground heat flux as the residual of the net radiation budget and turbulent fluxes.
"""
@kernel function compute_ground_heat_flux_kernel!(state, grid, ::AbstractSkinTemperature)
    i, j = @index(Global, NTuple)
    # compute flux terms
    R_net = state.surface_net_radiation[i, j]
    H_s = state.sensible_heat_flux[i, j]
    H_l = state.latent_heat_flux[i, j]
    # compute ground heat flux as residual of R_net and turbulent fluxes
    state.ground_heat_flux[i, j, 1] = R_net - H_s - H_l
end

"""
    update_skin_temperature_kernel!(state, grid, skinT::PrognosticSkinTemperature)

Diagnose the skin temperature implied by the current `ground_heat_flux` and `ground_temperature`.
"""
@kernel function update_skin_temperature_kernel!(state, grid, skinT::ImplicitSkinTemperature)
    i, j = @index(Global, NTuple)
    # get thickness of topmost soil/ground grid cell
    field_grid = get_field_grid(grid)
    Δz₁ = Δzᵃᵃᶜ(i, j, field_grid.Nz, field_grid)
    # get ground heat flux
    G = state.ground_heat_flux[i, j]
    # get ground temperature
    Tₛ = state.ground_temperature[i, j]
    # compute new skin temperature T₀ by setting G equal to the half-cell heat flux
    κₛ = skin_thermal_conductivity(i, j, grid, state, skinT)
    state.skin_temperature[i, j, 1] = Tₛ - G * Δz₁ / (2*κₛ)
end

# Kernel functions

@inline skin_thermal_conductivity(i, j, grid, state, skinT::ImplicitSkinTemperature) = skinT.κₛ

@inline skin_temperature(i, j, grid, state, ::AbstractSkinTemperature) = @inbounds state.skin_temperature[i, j]

@inline ground_temperature(i, j, grid, state, ::AbstractSkinTemperature) = @inbounds state.ground_temperature[i, j]

@inline ground_heat_flux(i, j, grid, state, ::AbstractSkinTemperature) = @inbounds state.ground_heat_flux[i, j]
