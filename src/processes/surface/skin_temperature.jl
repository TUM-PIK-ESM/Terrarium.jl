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
R_{\\text{net}}(T_s) + H_s(T_s) + H_l(T_s) - (1 - f_{\\text{snow}})\\, G(T_s, T_g) - f_{\\text{snow}}\\, S(T_s, T_{\\text{snow}}) = 0
```
where ``R_{\\text{net}}`` is the net radiation budget, ``H_s`` is the sensible heat flux, ``H_l`` is the latent
heat flux from sublimation and evapotranspiration, ``G`` is the conductive flux from the skin into the
snow-free ground (``T_g`` its temperature), ``S`` is the conductive flux from the skin into the top of the
snowpack over the snow-covered fraction (``T_{\\text{snow}}`` its temperature), and ``f_{\\text{snow}}``
is the snow-covered area fraction (``f_{\\text{snow}} = 0`` and ``S`` absent without snow). ``G`` and
``S`` are each computed from their own *unblended* conduction target (see
[`ground_thermal_interface`](@ref) and [`snow_conduction_interface`](@ref)) rather than a single
cover-fraction-weighted blend of ground and snow properties — this keeps the flux delivered to the
snow-covered fraction from being diluted (or, for small ``f_{\\text{snow}}``, from being replaced
outright) by the bare-ground share when computing the snowpack's own top-of-pack energy tendency (see
[`compute_skin_temperature_residual!`](@ref)). All fluxes follow the positive upwards convention: for
``G`` and ``S``, this means from the medium below towards the surface is positive; for the other fluxes,
this means from the surface towards the atmosphere is positive.

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
    return RootSolver(NF; max_iterations = 10)
end

"""
    $TYPEDSIGNATURES

Seed the prognostic `skin_temperature` with the current `ground_temperature` so the implicit
nonlinear solve starts from a physically sensible guess close to the root.
"""
function initialize!(state, grid, ::ImplicitSkinTemperature, args...)
    set!(state.skin_temperature, state.ground_temperature)
    return nothing
end

"""
    $TYPEDSIGNATURES

Return the conduction target `(Tg, κ, Δz)` for the snow-free ground: the uppermost ground layer's
`ground_temperature`, the assumed surface conductivity `κₛ`, and the top ground cell thickness. This is
always the *unblended* ground-only target, used both when there is no snow and (weighted by
`1 - f_snow`) for the bare-ground share of the skin-temperature solve when there is.
"""
@propagate_inbounds function ground_thermal_interface(i, j, grid, fields, skinT::ImplicitSkinTemperature)
    field_grid = get_field_grid(grid)
    Δz₁ = Δzᵃᵃᶜ(i, j, field_grid.Nz, field_grid)
    Tg = fields.ground_temperature[i, j]
    return (Tg, skinT.κₛ, Δz₁)
end

"""
    $TYPEDSIGNATURES

Return the conduction target `(Tsnow, κsnow, dsnow)` for the snow-covered fraction: the snow's own
(bulk) temperature, its thermal conductivity recovered from the density scheme, and its depth floored at
`min_conduction_thickness` (avoids a `0/0` conductance as the pack vanishes — the same floor already used
for the snow-base flux, see [`compute_snow_basal_heat_flux`](@ref)). Without snow (`snow === nothing`)
returns a finite placeholder (`κsnow = 0`) so the always-zero `f_snow` weight in
[`compute_skin_temperature_residual!`](@ref) eliminates this term without risking `0 * Inf = NaN`.
"""
@propagate_inbounds snow_conduction_interface(i, j, grid, fields, ::Nothing, constants::PhysicalConstants) =
    (zero(eltype(grid)), zero(eltype(grid)), one(eltype(grid)))
@propagate_inbounds function snow_conduction_interface(i, j, grid, fields, snow::AbstractSnow, constants::PhysicalConstants)
    ρ_snow = compute_snow_density(i, j, grid, fields, snow.density)
    κ_snow = compute_thermal_conductivity(snow, constants.material, ρ_snow)
    T_snow = snow_temperature(i, j, grid, fields, snow)
    d_snow = max(snow_depth(i, j, grid, fields, snow), snow.min_conduction_thickness)
    return (T_snow, κ_snow, d_snow)
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

Compute the ground heat flux from the surface net radiation and sensible/latent heat flux at grid cell
`i, j`: the flux the atmosphere-side physics *demands* of the surface, `G = R_net + H_s + H_l`. This is
the generic definition used directly by [`PrescribedSkinTemperature`](@ref) (which has no separate
conduction target to solve against, so the demanded and stored fluxes coincide by construction) and,
internally, as the atmosphere-side term of the [`ImplicitSkinTemperature`](@ref) residual (see
[`compute_skin_temperature_residual!`](@ref)) — but *not* as `ImplicitSkinTemperature`'s own
`ground_heat_flux` (see the more specific method below).
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

`ImplicitSkinTemperature`'s `ground_heat_flux` is the *explicit* snow-free-ground conductive flux
`G = 2κg(Tg − Ts)/Δzg`, evaluated directly from the current `skin_temperature` (whatever value the state
holds — a Newton iterate mid-solve, or the converged value) via [`ground_thermal_interface`](@ref). This
overrides the generic `R_net + H_s + H_l` method above: computing `G` explicitly from conduction (rather
than deriving it by algebraically inverting the atmosphere-side demand, as a previous implementation did)
keeps it an unambiguous, always-consistent per-bare-ground-area quantity at every iterate — not only in
the limit of exact convergence, and not conflated with the whole-cell atmosphere-side demand once a
separate snow-top conductive term also exists (see [`compute_skin_temperature_residual!`](@ref)).
"""
@propagate_inbounds function compute_ground_heat_flux(
        i, j, grid, fields,
        skinT::ImplicitSkinTemperature,
        ::AbstractSurfaceEnergyBalance
    )
    Tg, κg, Δzg = ground_thermal_interface(i, j, grid, fields, skinT)
    Ts = fields.skin_temperature[i, j]
    return 2 * κg * (Tg - Ts) / Δzg
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

Invert the (linear) area-weighted conduction relation for the implicit skin temperature `Ts` given the
atmosphere-side demanded flux `G` (`= R_net + H_s + H_l`), by equating `G` to the area-weighted sum of the
*unblended* ground and snow-top conductive fluxes,
`(1 − f_snow)·2κg(Tg − Ts)/Δzg + f_snow·2κsnow(Tsnow − Ts)/dsnow`. Since both conductive terms are affine
in `Ts` (only `G`, via `R_net`/`H`/`LE`, is nonlinear), this closed form is an *exact* solve of the
conduction side, not an approximation — the same linear-inversion structure used before the `G`/`S` split
(see [`compute_skin_temperature_residual!`](@ref)), just with the ground and snow terms kept unblended
rather than combined into a single cover-fraction-weighted conduction target. Without snow (`f_snow = 0`;
`κsnow`/`Tsnow` unused since their term is weighted out) this reduces exactly to `Ts = Tg − G·Δzg/(2κg)`.
"""
@inline function compute_skin_temperature(::ImplicitSkinTemperature, G, Tg, κg, Δzg, f_snow, Tsnow, κsnow, dsnow)
    Ag = 2 * κg / Δzg
    As = 2 * κsnow / dsnow
    A = (one(f_snow) - f_snow) * Ag + f_snow * As
    B = (one(f_snow) - f_snow) * Ag * Tg + f_snow * As * Tsnow
    Ts = (B - G) / A
    return Ts
end

"""
    $TYPEDSIGNATURES

Surface-energy-balance residual at grid cell `i, j`, in temperature space: `Ts_prev − Ts_implicit`, where
`Ts_implicit` is the exact conduction-side inverse (see [`compute_skin_temperature`](@ref)) of the
atmosphere-side demanded flux `G_demand = R_net(Ts_prev) + H(Ts_prev) + LE(Ts_prev)`. This keeps the
solver mechanics (and its tolerance semantics) unchanged from before the `G`/`S` split — the demanded flux
is recomputed here rather than read back from `ground_heat_flux`, because that field now stores the
*explicit* bare-ground conductive flux at the current `Ts` (see the `ImplicitSkinTemperature`-specific
[`compute_ground_heat_flux`](@ref)), not the atmosphere-side demand; the two coincide only at
convergence. `snow` also partitions the latent flux by area (see [`compute_surface_energy_fluxes!`](@ref)).
"""
@propagate_inbounds function compute_skin_temperature_residual!(
        out, i, j, grid, fields,
        skinT::ImplicitSkinTemperature,
        seb::AbstractSurfaceEnergyBalance,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        hydrology::Optional{AbstractSurfaceHydrology} = nothing,
        snow::Optional{AbstractSnow} = nothing
    )
    # Compute all fluxes at the current Ts (this includes the explicit ground-conduction flux, stored
    # into `ground_heat_flux`); `snow` partitions the latent flux by snow-covered fraction
    compute_surface_energy_fluxes!(out, i, j, grid, fields, seb, constants, atmos, hydrology, snow)
    R_net = out.surface_net_radiation[i, j, 1]
    H_s = out.sensible_heat_flux[i, j, 1]
    H_l = out.latent_heat_flux[i, j, 1]
    G_demand = compute_ground_heat_flux(skinT, R_net, H_s, H_l)
    Tg, κg, Δzg = ground_thermal_interface(i, j, grid, fields, skinT)
    f = snow_cover_fraction(i, j, grid, fields, snow)
    Tsnow, κsnow, dsnow = snow_conduction_interface(i, j, grid, fields, snow, constants)
    Ts_implicit = compute_skin_temperature(skinT, G_demand, Tg, κg, Δzg, f, Tsnow, κsnow, dsnow)
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
    Ts = solve!(out, (i, j), grid, fields, objective, skinT.solver, skinT, seb, args...)
    return Ts
end

# Kernels

@kernel function compute_ground_heat_flux_kernel!(out, grid, fields, skinT::AbstractSkinTemperature, args...)
    i, j = @index(Global, NTuple)

    # Compute and store ground heat flux
    out.ground_heat_flux[i, j, 1] = compute_ground_heat_flux(i, j, grid, fields, skinT, args...)
end
