"""
    $TYPEDEF

Represents skin temperature prescribed from input data.
"""
struct PrescribedSkinTemperature <: AbstractSkinTemperature end

variables(::PrescribedSkinTemperature) = (
    auxiliary(:ground_heat_flux, XY(), units=u"W/m^2", desc="Ground heat flux"),
    input(:skin_temperature, XY(), units=u"°C", desc="Longwave emission temperature of the land surface in °C"),
    input(:ground_temperature, XY(), units=u"°C", desc="Temperature of the uppermost ground or soil grid cell in °C")
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

variables(::PrognosticSkinTemperature) = (
    prognostic(:skin_temperature, XY(), units=u"°C", desc="Longwave emission temperature of the land surface in °C"),
    auxiliary(:ground_heat_flux, XY(), units=u"W/m^2", desc="Ground heat flux"),
    input(:ground_temperature, XY(), units=u"°C", desc="Temperature of the uppermost ground or soil grid cell in °C")
)

function compute_auxiiliary!(state, model, skinT::PrognosticSkinTemperature)
    grid = get_grid(model)
    rad = get_radiative_fluxes(model)
    tur = get_turbulent_fluxes(model)
    launch!(grid, compute_ground_heat_flux!, state, grid, skinT, rad, tur)
end

function compute_tendencies!(state, model, ::PrognosticSkinTemperature)
    grid = get_grid(model)
    launch!(grid, compute_skin_temperature_tendency!, state, grid, skinT)
end

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

@kernel function compute_skin_temperature_tendency!(
    state, grid,
    skinT::PrognosticSkinTemperature
)
    i, j = idx = @index(Global, NTuple)
    # compute and accumulate tendency
    state.tendencies.skin_temperature[i, j] += skin_temperature_tendency(idx, state, grid, skinT)
end

# Kernel functions

@inline function skin_temperature_tendency(idx, state, grid, skinT::PrognosticSkinTemperature)
    i, j = idx
    # get thickness of topmost soil/ground grid cell
    field_grid = get_field_grid(grid)
    Δz₁ = Δzᵃᵃᶜ(i, j, field_grid.Nz, field_grid)
    # get ground heat flux
    G = state.ground_heat_flux[i, j]
    # compute correction ΔTₛ by setting G equal to the half-cell heat flux;
    # note that this is technically not a tendency and should *not* be multiplied by Δt!
    κₛ = skin_thermal_conductivity(idx, state, skinT)
    ΔTₛ = -G * Δz₁ / (2*κₛ)
    return ΔTₛ
end

@inline skin_thermal_conductivity(idx, state, skinT::PrognosticSkinTemperature) = skinT.κₛ

@inline @inbounds skin_temperature(idx, state, ::AbstractSkinTemperature) = state.skin_temperature[idx]

@inline @inbounds ground_heat_flux(idx, state, ::AbstractSkinTemperature) = state.ground_heat_flux[idx]
