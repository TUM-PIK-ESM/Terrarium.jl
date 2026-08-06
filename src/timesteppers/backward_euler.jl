"""
    $TYPEDEF

Backward-Euler implicit time stepping scheme using ClimaTimeSteppers.jl `NewtonsMethod`
for the nonlinear solve and Enzyme forward-mode for Jacobian computation.

The implicit equation solved at each step is:

    G(u) = u − u_prev − Δt · f(u) = 0

where `f` is the tendency map (closure! + compute_tendencies!).

Properties:
$(TYPEDFIELDS)
"""
struct BackwardEuler{NF} <: AbstractTimeStepper{NF}
    "CTS Newton solver configuration"
    newton_method::ClimaTimeSteppers.NewtonsMethod

    "Timestep size in seconds"
    Δt::NF
end

"""
    BackwardEuler(::Type{NF}; Δt, max_iters, kwargs...)

Construct a `BackwardEuler` timestepper for numeric type `NF`.

# Arguments
- `max_iters`: maximum Newton iterations (default: 5)
- `Δt`: timestep size in seconds (default: 300.0)
"""
function BackwardEuler(
        ::Type{NF};
        Δt::Real = 300.0,
        max_iters::Int = 5,
        kwargs...
    ) where {NF}
    newton_method = ClimaTimeSteppers.NewtonsMethod(;
        max_iters,
        update_j = ClimaTimeSteppers.UpdateEvery(ClimaTimeSteppers.NewNewtonIteration),
        kwargs...,
    )
    return BackwardEuler(newton_method, convert(NF, Δt))
end

timestepping(::BackwardEuler) = Implicit()

default_dt(be::BackwardEuler) = be.Δt

is_adaptive(be::BackwardEuler) = false

"""
    $TYPEDSIGNATURES

Write the flat state vector `u` into the interior of the prognostic field.
"""
@inline function u2state!(state, u)
    field = state.prognostic.internal_energy
    interior(field) .= reshape(u, size(interior(field)))
    return nothing
end

"""
    $TYPEDSIGNATURES

Read the interior of the prognostic field as a flat vector.
"""
@inline function state2u(state)
    return vec(interior(state.prognostic.internal_energy))
end

"""
    $TYPEDSIGNATURES

Compute the tendency map `f(u)` — the time derivative of internal energy — at the
current state. Calls `fill_halo_regions!` to ensure boundary conditions track the
current interior temperature.
"""
@inline function compute_f!(state, grid, soil, constants)
    fill_halo_regions!(state)
    Terrarium.reset_tendencies!(state)
    Terrarium.closure!(state, grid, soil, constants)
    Terrarium.compute_tendencies!(state, grid, soil, constants)
    return nothing
end

# CTS convention: f_CTS(x) = x_prev + Δt·f(x) − x  (note the sign flip)
# This gives J_CTS = Δt·J_f − I, and Newton update x ← x − J_CTS⁻¹·f_CTS
#
# The residual is computed through Fields so that Enzyme can trace the full
# U → T → tendencies → residual chain. The flat-vector conversion happens
# only at the CTS interface.
@inline function backward_euler_residual!(res, u, u_prev, state, grid, soil, constants, Δt)
    u2state!(state, u)
    compute_f!(state, grid, soil, constants)
    # Compute residual as a Field operation, then extract to flat vector for CTS
    tendencies_field = state.tendencies.internal_energy
    u_field = state.prognostic.internal_energy
    res .= u_prev .+ Δt .* vec(interior(tendencies_field)) .- vec(interior(u_field))
    return nothing
end

"""
    $TYPEDSIGNATURES

Compute the Jacobian `W = Δt · J_f − I` using Enzyme forward-mode automatic
differentiation. Each column of J_f is computed as a Jacobian-vector product (JVP)
via `Enzyme.autodiff(Forward, ...)`, following the same pattern as
`examples/autodiff/differentiating_terrarium.jl`.

The JVP computation goes through Fields via `Duplicated(state, dstate)`, which is
the same pattern that works for sensitivity analysis in the rest of Terrarium.
"""
@inline function compute_jacobian!(J, u, u_prev, dstate, state, grid, soil, constants, Δt)
    N_z = length(u)
    u2state!(state, u)
    compute_f!(state, grid, soil, constants)
    J_f = zeros(eltype(u), N_z, N_z)
    shadow_u = interior(dstate.prognostic.internal_energy)
    for col in 1:N_z
        shadow_u .= zero(eltype(u))
        shadow_u[1, 1, col] = one(eltype(u))
        Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.Forward),
            compute_f!, Enzyme.Const,
            Enzyme.Duplicated(state, dstate),
            Enzyme.Const(grid), Enzyme.Const(soil), Enzyme.Const(constants),
        )
        J_f[:, col] .= vec(interior(dstate.tendencies.internal_energy))
    end
    # Form W = Δt · J_f − I  (CTS convention)
    J .= Δt .* J_f
    for i in 1:N_z
        J[i, i] -= one(eltype(u))
    end
    return nothing
end

function timestep!(
        integrator::ModelIntegrator,
        timestepper::BackwardEuler,
        Δt,
        names::Tuple,
    )
    state = integrator.state
    grid = get_grid(integrator.model)
    soil = integrator.model.soil
    constants = integrator.model.constants

    # Extract state as flat vector (copy, not view — CTS mutates x in-place)
    u = copy(state2u(state))
    u_prev = copy(u)

    # Allocate Enzyme shadow state for JVP seeding
    dstate = Enzyme.make_zero(state)

    # Build residual and Jacobian closures
    f! = (res, x) -> backward_euler_residual!(res, x, u_prev, state, grid, soil, constants, Δt)
    j! = (J, x) -> compute_jacobian!(J, x, u_prev, dstate, state, grid, soil, constants, Δt)

    # Allocate Jacobian prototype (dense for now; CTS handles LU factorization)
    N_z = length(u)
    j_prototype = zeros(eltype(u), N_z, N_z)

    # Allocate CTS Newton solver cache
    cache = ClimaTimeSteppers.allocate_cache(timestepper.newton_method, u, j_prototype)

    # Solve the nonlinear system
    ClimaTimeSteppers.solve_newton!(timestepper.newton_method, cache, u, f!, j!)

    # Write solution back into field
    u2state!(state, u)

    return nothing
end
