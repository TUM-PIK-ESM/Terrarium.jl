# Prescribed skin temperature

"""
    $TYPEDEF

Simple scheme for prescribed skin temperatures from input variables.

Properties:
$FIELDS
"""
@parameterized @kwdef struct PrescribedSkinTemperature{NF} <: AbstractSkinTemperature{NF}
    "Assumed thermal conductivity at the surface"
    @param κₛ::NF = 1.0 (units = u"W/m/K", bounds = Positive)
end

PrescribedSkinTemperature(::Type{NF}; kwargs...) where {NF} = PrescribedSkinTemperature{NF}(; kwargs...)

## Top-level interface methods

variables(::PrescribedSkinTemperature) = (
    auxiliary(:ground_heat_flux, XY(), units = u"W/m^2", desc = "Ground heat flux"),
    input(:skin_temperature, XY(), units = u"°C", desc = "Longwave emission temperature of the land surface in °C"),
)

@inline compute_auxiliary!(state, grid, ::PrescribedSkinTemperature, args...) = nothing

@inline compute_tendencies!(state, grid, ::PrescribedSkinTemperature, args...) = nothing

## Kernel functions

@propagate_inbounds compute_skin_temperature(i, j, grid, fields, skinT::PrescribedSkinTemperature) = fields.skin_temperature[i, j]

# The skin temperature is prescribed (an input field), so there is no implicit solve: the fused
# surface-energy-balance kernel just evaluates the fluxes from the prescribed skin temperature.
@propagate_inbounds solve_skin_temperature!(out, i, j, grid, fields, ::PrescribedSkinTemperature, seb, snow, seb_args...) = nothing

# Implicit skin temperature

"""
    $TYPEDEF

Scheme for an implicit skin temperature ``T_s`` satisfying:
```math
R_{\\text{net}}(T_s) + H_s(T_s) + H_l(T_s) - G(T_s, T_g) = 0 
```
where ``R_{\\text{net}}`` is the net radiation budget, ``H_s`` is the sensible heat flux, ``H_l`` is the latent
heat flux from sublimation and evapotranspiration, ``G`` is the ground heat flux, and ``T_g`` is the ground
temperature, or temperature of the uppermost subsurface (soil or snow) layer. All fluxes follow the positive upwards 
convention. For ``G``, this means from the deeper soil towards the surface is positive. For the other fluxes, this 
means from the surface towards the atmosphere is positive. 

Properties:
$FIELDS
"""
@parameterized struct ImplicitSkinTemperature{NF, Solver} <: AbstractSkinTemperature{NF}
    "Assumed thermal conductivity at the surface"
    @param κₛ::NF (units = u"W/m/K", bounds = Positive)

    "Numerical solver for the implicit skin temperature"
    solver::Solver
end

ImplicitSkinTemperature(::Type{NF}; κₛ::NF = NF(1.0), solver = default_skin_temperature_solver(NF)) where {NF} = ImplicitSkinTemperature{NF, typeof(solver)}(κₛ, solver)

"""
    $TYPEDSIGNATURES

Construct the default solver for the implicit skin temperature: a Newton root-finder
([`RootSolver`](@ref) backed by RootSolvers.jl) with a small iteration budget.
"""
function default_skin_temperature_solver(::Type{NF}) where {NF}
    # Default to a Newton root-finder (via RootSolvers.jl) with a small iteration budget
    return RootSolver(NF; max_iterations = 5)
end

"""
    $TYPEDSIGNATURES

Seed the prognostic `skin_temperature` with the current `ground_temperature` so the implicit
nonlinear solve starts from a physically sensible guess close to the root.
"""
function initialize!(state, grid, ::ImplicitSkinTemperature, args...)
    interior(state.skin_temperature) .= interior(state.ground_temperature)   # Now set! for better Reactant compat, set! might make a CPU detour
    return nothing
end

"""
    $TYPEDSIGNATURES

Compute the implicit update of the skin temperature from the sub-surface temperature `Tg`, ground heat
flux `G`, half-cell distance `Δz`, and effective conductivity `κ`, by equating `G` to the half-cell
conductive flux `2κ(Tg − Ts)/Δz`.
"""
@inline function compute_skin_temperature(::ImplicitSkinTemperature, Tg, G, Δz, κ)
    Ts = Tg - G * Δz / (2 * κ)
    return Ts
end

"""
    $TYPEDSIGNATURES

Return the effective conduction target `(Tg, κ, Δz)` seen by the skin-temperature solve: the temperature
of the medium directly beneath the skin, its thermal conductivity, and the half-cell thickness. Without
snow (`snow === nothing`) these are the uppermost ground layer's `ground_temperature`, the assumed
surface conductivity `κₛ`, and the top ground cell thickness — i.e. the snow-free behavior is unchanged.
"""
@propagate_inbounds function ground_thermal_interface(i, j, grid, fields, skinT::ImplicitSkinTemperature, snow::Nothing, constants = nothing)
    field_grid = get_field_grid(grid)
    Δz₁ = Δzᵃᵃᶜ(i, j, field_grid.Nz, field_grid)
    Tg = fields.ground_temperature[i, j]
    return (Tg, skinT.κₛ, Δz₁)
end

"""
    $TYPEDSIGNATURES

Snow-covered variant: the conduction target is the snow-cover-fraction-weighted blend of the snow layer
and the underlying ground, `x_eff = (1 − f_snow)·x_ground + f_snow·x_snow` for the temperature,
conductivity, and half-cell thickness.
"""
@propagate_inbounds function ground_thermal_interface(i, j, grid, fields, skinT::ImplicitSkinTemperature, snow::AbstractSnow, constants::PhysicalConstants)
    field_grid = get_field_grid(grid)
    Δz_ground = Δzᵃᵃᶜ(i, j, field_grid.Nz, field_grid)
    f = snow_cover_fraction(i, j, grid, fields, snow)
    # bulk snow thermal conductivity recovered lazily from the density scheme (not stored as a field)
    ρ_snow = compute_snow_density(i, j, grid, fields, snow.density)
    κ_snow = compute_thermal_conductivity(snow, constants.material, ρ_snow)
    Tg = (one(f) - f) * fields.ground_temperature[i, j] + f * snow_temperature(i, j, grid, fields, snow)
    κ = (one(f) - f) * skinT.κₛ + f * κ_snow
    Δz = (one(f) - f) * Δz_ground + f * snow_depth(i, j, grid, fields, snow)
    return (Tg, κ, Δz)
end

"""
    $TYPEDSIGNATURES

Compute the ground heat flux `G` that closes the surface energy balance. With all fluxes
positive upward (aligned with `+z`), the energy arriving at the skin from below must balance
the radiative and turbulent losses above, so `G = R_net + H_s + H_l`.
"""
@inline function compute_ground_heat_flux(::AbstractSkinTemperature, R_net, H_s, H_l)
    G = R_net + H_s + H_l
    return G
end

## Top-level interface methods

variables(::ImplicitSkinTemperature) = (
    prognostic(:skin_temperature, XY(), units = u"°C", desc = "Longwave emission temperature of the land surface in °C"),
    auxiliary(:ground_heat_flux, XY(), units = u"W/m^2", desc = "Ground heat flux"),
    input(:ground_temperature, XY(), units = u"°C", desc = "Temperature of the uppermost ground or soil grid cell in °C"),
)

""" $TYPEDSIGNATURES """
@inline function compute_auxiliary!(
        state, grid,
        skinT::ImplicitSkinTemperature,
        seb::AbstractSurfaceEnergyBalance,
        args...
    )
    compute_ground_heat_flux!(state, grid, skinT, seb)
    return nothing
end

"""
    $TYPEDSIGNATURES

Compute `skin_temperature` from the current `ground_heat_flux` using the (optionally snow-aware)
conduction target. `constants` supplies the material properties for the conduction blend; passing a
`snow` component (default `nothing`) selects the snow-covered [`ground_thermal_interface`](@ref), while
`snow === nothing` recovers the snow-free soil conduction.
"""
function compute_skin_temperature!(
        state, grid,
        skinT::ImplicitSkinTemperature,
        constants::PhysicalConstants,
        snow::Optional{AbstractSnow} = nothing,
    )
    out = prognostic_fields(state, skinT)
    # merge the (optional) snow thermal auxiliaries so the snow-aware conduction target can be evaluated
    fields = merge(get_fields(state, skinT; except = out), get_fields(state, snow))
    launch!(grid, XY, compute_skin_temperature_kernel!, out, fields, skinT, constants, snow)
    return nothing
end

"""
    $TYPEDSIGNATURES

Compute `ground_heat_flux` as the residual ``R_\\text{net} + H_s + H_l`` (all fluxes positive upward).
"""
function compute_ground_heat_flux!(
        state, grid,
        skinT::AbstractSkinTemperature,
        seb::AbstractSurfaceEnergyBalance
    )
    out = auxiliary_fields(state, skinT)
    fields = get_fields(state, seb; except = out)
    launch!(grid, XY, compute_ground_heat_flux_kernel!, out, fields, skinT, seb)
    return nothing
end

## Kernel functions

"""
    $TYPEDSIGNATURES

Estimate the skin temperature from the current `ground_heat_flux` at grid cell `i, j`, using the
effective conduction target of the medium below the skin. The optional `snow` argument enables the
snow-covered blend (see [`ground_thermal_interface`](@ref)); with `snow === nothing` this is the snow-free soil
conduction.
"""
@propagate_inbounds function compute_skin_temperature(
        i, j, grid, fields,
        skinT::ImplicitSkinTemperature,
        constants::PhysicalConstants,
        snow::Optional{AbstractSnow} = nothing
    )
    Tg, κ, Δz = ground_thermal_interface(i, j, grid, fields, skinT, snow, constants)
    G = fields.ground_heat_flux[i, j]
    return compute_skin_temperature(skinT, Tg, G, Δz, κ)
end

"""
    $TYPEDSIGNATURES

Compute the ground heat flux from the surface net radiation and sensible/latent heat flux at grid cell `i, j`.
"""
@propagate_inbounds function compute_ground_heat_flux(
        i, j, grid, fields,
        skinT::AbstractSkinTemperature,
        ::AbstractSurfaceEnergyBalance
    )
    # Get individual flux terms
    R_net = fields.surface_net_radiation[i, j]
    H_s = fields.sensible_heat_flux[i, j]
    H_l = fields.latent_heat_flux[i, j]
    # Compute ground heat flux
    G = compute_ground_heat_flux(skinT, R_net, H_s, H_l)
    return G
end

"""
    $TYPEDSIGNATURES

Per-cell mutating variant used by the fused surface-energy-balance kernel: store the ground heat flux
into the auxiliary output field `out`.
"""
@propagate_inbounds function compute_ground_heat_flux!(out, i, j, grid, fields, skinT::AbstractSkinTemperature, seb::AbstractSurfaceEnergyBalance)
    out.ground_heat_flux[i, j, 1] = compute_ground_heat_flux(i, j, grid, fields, skinT, seb)
    return nothing
end

"""
    $TYPEDSIGNATURES

Recompute surface energy fluxes using [`compute_surface_energy_fluxes!`](@ref) and update the skin temperature
based on the resulting ground heat flux.
"""
@propagate_inbounds function compute_skin_temperature!(
        out, i, j, grid, fields,
        skinT::ImplicitSkinTemperature,
        seb::AbstractSurfaceEnergyBalance,
        constants::PhysicalConstants,
        snow::Optional{AbstractSnow} = nothing,
        seb_args...
    )
    # Compute all fluxes based on current skin temperature (this includes ground heat flux already);
    # `snow` (after evtr = nothing) partitions the latent flux by snow-covered fraction
    compute_surface_energy_fluxes!(out, i, j, grid, fields, seb, constants, seb_args..., nothing, snow)
    out.skin_temperature[i, j, 1] = compute_skin_temperature(i, j, grid, fields, skinT, constants, snow)
    return nothing
end

"""
    $TYPEDSIGNATURES

Same as [`compute_skin_temperature!`](@ref) but returns the residual instead of updating the `skin_temperature` `Field`.
"""
@propagate_inbounds function compute_skin_temperature_residual!(
        out, i, j, grid, fields,
        skinT::ImplicitSkinTemperature,
        seb::AbstractSurfaceEnergyBalance,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        evtr::Optional{AbstractEvapotranspiration} = nothing,
        snow::Optional{AbstractSnow} = nothing
    )
    # Compute all fluxes based on current skin temperature (this includes ground heat flux already);
    # `snow` (after evtr = nothing) partitions the latent flux by snow-covered fraction
    compute_surface_energy_fluxes!(out, i, j, grid, fields, seb, constants, atmos, evtr, snow)
    # Compute skin temperature using the (optionally snow-aware) conduction target
    Ts_implicit = compute_skin_temperature(i, j, grid, fields, skinT, constants, snow)
    Ts_prev = out.skin_temperature[i, j, 1]
    return Ts_prev - Ts_implicit
end

"""
    $TYPEDSIGNATURES

Run a full nonlinear solve to determine the `skin_temperature` at grid cell `i, j` that solves the surface energy balance.
"""
@propagate_inbounds function solve_skin_temperature!(
        out, i, j, grid, fields,
        skinT::ImplicitSkinTemperature,
        seb::AbstractSurfaceEnergyBalance,
        args...
    )
    objective = ObjectiveFunction(compute_skin_temperature_residual!, :skin_temperature)
    # `snow` is forwarded (after `seb`) so the residual's conduction target can be snow-aware
    return solve!(out, (i, j), grid, fields, objective, skinT.solver, skinT, seb, args...)
end

# Kernels

@kernel function compute_skin_temperature_kernel!(out, grid, fields, skinT::ImplicitSkinTemperature, args...)
    i, j = @index(Global, NTuple)

    out.skin_temperature[i, j, 1] = compute_skin_temperature(i, j, grid, fields, skinT, args...)
end

@kernel function compute_ground_heat_flux_kernel!(out, grid, fields, skinT::AbstractSkinTemperature, args...)
    i, j = @index(Global, NTuple)

    # Compute and store ground heat flux
    out.ground_heat_flux[i, j, 1] = compute_ground_heat_flux(i, j, grid, fields, skinT, args...)
end
