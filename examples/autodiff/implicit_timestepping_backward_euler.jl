# # Implicit timestepping with ClimaTimeSteppers.jl Newton + Enzyme Jacobian
#
# Phase 3 of the implicit-timestepping experiment. A backward-Euler step requires solving
#
#   G(u) = u − u_prev − Δt · f(u) = 0
#
# where `f` is the Terrarium tendency map. This example uses [ClimaTimeSteppers.jl](https://github.com/CliMA/ClimaTimeSteppers.jl)
# `NewtonsMethod` as the nonlinear solver and [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl)
# forward-mode for Jacobian computation.
#
# The key design decisions are:
#
#   1. **Residual through Fields**: The residual is computed via Oceananigans `Field` objects
#      (not flat vectors) so that Enzyme can trace the full U→T→tendencies→residual chain.
#      Flat-vector conversion happens only at the CTS interface.
#
#   2. **Jacobian via column-by-column JVPs**: `Enzyme.jacobian(Forward, ...)` fails with
#      `EnzymeMutabilityException` because the closure captures mutable state. The working pattern
#      is `Enzyme.autodiff(Forward, compute_f!, Const, Duplicated(state, dstate), ...)` with
#      one-hot seeding per column — the same pattern used in `differentiating_terrarium.jl`.
#
#   3. **CTS used minimally**: Only `NewtonsMethod`, `allocate_cache`, `solve_newton!` from
#      ClimaTimeSteppers. No `ODEFunction`, `ClimaODEFunction`, or `ODEProblem`.
#
#   4. **`vec(interior(field))` not `parent(field)`**: Excludes halos; uses Oceananigans-recommended
#      access pattern. `copy()` is required because CTS mutates `x` in-place via `x .-= Δx`.
#
#   5. **`I` (UniformScaling) cannot be broadcast**: Use explicit `J[i,i] -= one(eltype(u))` loop
#      instead of `J .= Δt .* J_f .- I`.
#
# All helper functions (`compute_f!`, `u2state!`, `state2u`, `backward_euler_residual!`,
# `compute_jacobian!`) are defined in `src/timesteppers/backward_euler.jl` and used here
# via qualified names.

using Terrarium
import ClimaTimeSteppers
import Oceananigans
using Enzyme
using LinearAlgebra

# ## Model setup
#
# Single column, uniform spacing, prescribed surface temperature of 1°C at the top boundary.

FT = Float64
grid = ColumnGrid(CPU(), FT, UniformSpacing(), 2)
initializer = SoilInitializer(eltype(grid))
model = SoilModel(grid; initializer, timestepper = BackwardEuler(FT; Δt = 300.0))
bcs = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator = initialize(model; boundary_conditions = bcs)
Terrarium.initialize!(integrator) # reset to t = 0
state = integrator.state
N_z = size(interior(state.prognostic.internal_energy), 3)

println("Grid Nz: $N_z")
println("Initial energy range: $(extrema(vec(interior(state.prognostic.internal_energy))))")

# ## Key building blocks (from `src/timesteppers/backward_euler.jl`)
#
# The following helper functions are defined in the source and used below:
#
# - `Terrarium.compute_f!(state, grid, soil, constants)`: combined tendency map
#   (fill_halo_regions! + reset_tendencies! + closure! + compute_tendencies!).
# - `Terrarium.u2state!(state, u)`: write flat vector `u` into interior of prognostic field.
# - `Terrarium.state2u(state)`: read interior of prognostic field as flat vector.
# - `Terrarium.backward_euler_residual!(res, u, u_prev, ...)`: CTS-compatible residual
#   `f_CTS(x) = x_prev + Δt·f(x) − x`.
# - `Terrarium.compute_jacobian!(J, u, u_prev, ...)`: Jacobian `W = Δt·J_f − I` via
#   column-by-column Enzyme forward-mode JVPs.

# ## Single backward-Euler timestep via `timestep!`
#
# The `BackwardEuler` timestepper's `timestep!` method (defined in `backward_euler.jl`)
# orchestrates the full solve: copy state → flat vector, allocate Enzyme shadow,
# build closures, call `CTS.solve_newton!`, write solution back.

println("\n=== Backward Euler timestep ===")
u_before = copy(vec(interior(state.prognostic.internal_energy)))
println("Energy before: $(extrema(u_before))")

@time timestep!(integrator)

u_after = copy(vec(interior(state.prognostic.internal_energy)))
println("Energy after:  $(extrema(u_after))")
println("Energy change: $(maximum(abs.(u_after .- u_before)))")

# ## Forward-Euler reference for comparison
#
# Re-initialize and run one forward-Euler step with the same Δt for comparison.

grid_fe = ColumnGrid(CPU(), FT, UniformSpacing(), 1)
model_fe = SoilModel(grid_fe; initializer = SoilInitializer(eltype(grid_fe)))
bcs_fe = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator_fe = initialize(model_fe; boundary_conditions = bcs_fe)
Terrarium.initialize!(integrator_fe)
Δt_fe = integrator_fe.model.timestepper.Δt

println("\n=== Forward Euler reference (Δt = $Δt_fe) ===")
println("Energy before: $(extrema(vec(interior(integrator_fe.state.prognostic.internal_energy))))")
@time timestep!(integrator_fe)
println("Energy after:  $(extrema(vec(interior(integrator_fe.state.prognostic.internal_energy))))")

# ## Multi-step comparison
#
# Run both methods for several timesteps and compare.

println("\n=== Multi-step comparison ===")
Terrarium.initialize!(integrator)    # reset backward Euler
Terrarium.initialize!(integrator_fe) # reset forward Euler

nsteps = 10
for i in 1:nsteps
    timestep!(integrator)
    timestep!(integrator_fe)
end

T_be = vec(interior(integrator.state.temperature))
T_fe = vec(interior(integrator_fe.state.temperature))
println("Backward Euler T range: $(extrema(T_be))")
println("Forward Euler  T range: $(extrema(T_fe))")
println("Max |ΔT| = $(maximum(abs.(T_be .- T_fe)))")
println("Relative error = $(norm(T_be .- T_fe) / norm(T_fe))")
