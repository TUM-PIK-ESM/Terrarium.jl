# # Implicit timestepping with NonlinearSolve.jl
#
# Experiment: solve the backward-Euler step of the Terrarium soil heat equation with
# [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl) instead of the Ariadne.jl
# prototype used in `differentiating_terrarium.jl`. A backward-Euler step requires solving
#
#   G(u) = u − u_prev − Δt·f(u) = 0
#
# where f is the Terrarium tendency map. The Jacobian of G is W = I − Δt·J_f with J_f the
# (tridiagonal) Jacobian of the heat-equation tendency w.r.t. internal energy.
#
# The key difficulty is that the default AD backend of NonlinearSolve (ForwardDiff) cannot
# differentiate through the Terrarium state wrapper: the state arrays are Float64, so dual
# numbers (and sparsity tracers) cannot be written into them. We therefore compare solver
# configurations that avoid ForwardDiff:
#
#   1. Manual Enzyme JVP (`jvp` callback) + Krylov GMRES, with and without Eisenstat–Walker forcing
#   2. `jvp_autodiff = AutoEnzyme` + GMRES (experimental: DI ↔ Enzyme in-place path)
#   3. Manual colored Enzyme Jacobian (3 multi-hot JVPs) + CSC tridiagonal prototype + KLU
#   4a. Automatic sparsity detection Route A: `DenseSparsityDetector(AutoEnzyme)` + KLU
#       (experimental; needs SparseMatrixColorings for the colored sparse Jacobian)
#   4b. Fully automatic sparsity exploitation with a working backend:
#       `DenseSparsityDetector(AutoFiniteDiff)` + colored FD Jacobian + KLU
#   4c. Route A with Enzyme done manually: dense Enzyme Jacobian → sparse pattern
#       (detected once) + colored FD Jacobian + KLU
#   5. Automatic sparsity detection Route B: `TracerLocalSparsityDetector` pattern from a
#       tracer-typed state + colored FD Jacobian + KLU
#   6. `AutoFiniteDiff` dense Jacobian + LU (fallback baseline)
#   7. Default `solve(prob)` polyalgorithm (its default ForwardDiff backend cannot handle
#      the Float64 state wrapper, but the polyalgorithm recovers via fallback methods)
#
# All residual formulations are built exclusively from Terrarium's own methods
# (`reset_tendencies!` → `closure!` → `compute_tendencies!`); no physics is re-implemented.
# For the tracer-based detection (variant 5), the *container* is re-typed instead of the
# physics: a tracer-typed copy of the `StateVariables` is constructed with Oceananigans'
# `Field{LX, LY, LZ}(grid, T)` constructor, and the unmodified `compute_f!` runs on it
# (tracing the Float64 state directly fails at the `u → state` write with a `TypeError`).
#
# Correctness is checked against the Ariadne.jl Newton–Krylov reference implementation.
#
# Note on the DI↔Enzyme failures (variants 2, 4a), informed by the Enzyme FAQ
# (https://enzymead.github.io/Enzyme.jl/dev/faq/): passing a runtime-activity-enabled mode
# through DI — `AutoEnzyme(mode = Enzyme.set_runtime_activity(Enzyme.Forward))` — does
# eliminate the `EnzymeRuntimeActivityError` in the pushforward path, but the resulting
# JVP is *silently wrong*: DI marks the state-holding parameter `p` as `Constant`, while
# the residual's active dataflow runs through `p.state` — precisely the FAQ's "activity of
# temporary storage" trap (values inside a `Const` argument are definitionally constant,
# so the derivative loses the state path). The dense-Jacobian path additionally crashes
# inside Enzyme. Hence manual callbacks with explicit `Duplicated(state, dstate)` (the
# FAQ's recommended fix) remain the only correct Enzyme route; variant 4c uses them to
# provide Enzyme-based sparsity detection.

using Terrarium
using Enzyme
using Enzyme: Forward, set_runtime_activity
using ADTypes
using DifferentiationInterface: DenseSparsityDetector
using Ariadne
using BenchmarkTools
using CUDA
using LinearAlgebra
using LinearSolve
using Krylov
using NonlinearSolve
import Oceananigans
using SparseArrays
using SparseConnectivityTracer
using SparseMatrixColorings

# ## Model setup
#
# Identical to the Ariadne prototype: single column, uniform spacing, prescribed surface
# temperature of 1°C at the top boundary.

FT = Float64
grid = ColumnGrid(CPU(), FT, UniformSpacing(), 1)
initializer = SoilInitializer(eltype(grid))
model = SoilModel(grid; initializer)
bcs = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator = initialize(model, boundary_conditions = bcs)
Terrarium.initialize!(integrator) # reset to t = 0
state = integrator.state
Δt = integrator.model.timestepper.Δt
N_z = size(interior(state.internal_energy), 3)

# Combined closure + tendency map so the full U → T → ∂U/∂t chain is visible to Enzyme.
# Note that, as in the reference example, `compute_f!` does not re-fill halo regions:
# the top ghost cell of `temperature` keeps the value set by `fill_halo_regions!` during
# initialization (2·T_bc − T_init(N_z) for the top ValueBoundaryCondition) and all other
# ghost cells keep their initial (zero) values. These frozen-ghost semantics are part of
# the residual operator and are exactly reproduced by the column-local residual below.
function compute_f!(state, grid, soil, constants)
    Terrarium.reset_tendencies!(state)
    Terrarium.closure!(state, grid, soil, constants)
    Terrarium.compute_tendencies!(state, grid, soil, constants)
    return nothing
end

u_prev = copy(vec(interior(state.prognostic.internal_energy)))
u0 = copy(u_prev)

# Solver tolerances: internal energy is O(10⁷) J/m³, so an absolute tolerance of 1e-2
# on the residual is physically tight. reltol = 0 forces pure absolute termination.
abstol = FT(1.0e-2)
reltol = FT(0)

# ## Vector ↔ state bridge and residual
#
# Shared parameter bundle for all residual/Jacobian callbacks. `dstate` is the preallocated
# Enzyme shadow of `state`; see the seeding protocol in `jvp_nls!`.

p = (
    state = state,
    dstate = make_zero(state),
    grid = model.grid,
    soil = model.soil,
    constants = model.constants,
    u_prev = u_prev,
    Δt = Δt,
)

u2state!(state, u) = interior(state.prognostic.internal_energy) .= reshape(u, size(interior(state.prognostic.internal_energy)))

function G_nls!(res, u, p)
    u2state!(p.state, u)
    compute_f!(p.state, p.grid, p.soil, p.constants)
    res .= u .- p.u_prev .- p.Δt .* vec(interior(p.state.tendencies.internal_energy))
    return nothing
end

# ## Manual Enzyme Jacobian-vector product (variants 1)
#
# J_G·v = v − Δt·J_f·v. The forward-mode seeding protocol: zero the shadow of the
# prognostic variable and seed it with v. All other shadow arrays are handled correctly
# without manual zeroing:
#   - `reset_tendencies!` zeroes the (shadow) tendencies as the derivative of `x .= 0`,
#   - `closure!` *overwrites* the shadow closure variables (temperature, liquid fraction),
#   - ghost-cell shadows remain zero from the initial `make_zero`, which is correct since
#     the ghost values are constant w.r.t. u.
function jvp_nls!(Jv, v, u, p)
    u2state!(p.state, u)
    shadow_u = interior(p.dstate.prognostic.internal_energy)
    shadow_u .= 0
    shadow_u .= reshape(v, size(shadow_u))
    Enzyme.autodiff(
        set_runtime_activity(Forward), compute_f!, Const,
        Duplicated(p.state, p.dstate),
        Const(p.grid), Const(p.soil), Const(p.constants),
    )
    Jv .= v .- p.Δt .* vec(interior(p.dstate.tendencies.internal_energy))
    return nothing
end

# ## Manual colored Enzyme Jacobian (variant 3)
#
# J_f is exactly tridiagonal: the freeze-curve closure is pointwise (diagonal) and heat
# conduction couples nearest neighbors only (κ at a face depends on the liquid fraction
# of the two adjacent cells). Columns of the same color (color(j) = mod(j-1, 3)) are
# structurally orthogonal, so J_f is recovered exactly from 3 JVPs with multi-hot seeds:
# for row i, exactly one of the band columns {i−1, i, i+1} has color c, hence
# (J·v_c)[i] = J[i, j] for that column j.

tridiag_color(j) = mod(j - 1, 3) + 1

function jac_nls!(JW, u, p)
    u2state!(p.state, u)
    nzv = nonzeros(JW)
    shadow_u = interior(p.dstate.prognostic.internal_energy)
    for c in 1:3
        shadow_u .= 0
        for j in 1:N_z
            shadow_u[1, 1, j] = ifelse(tridiag_color(j) == c, one(FT), zero(FT))
        end
        Enzyme.autodiff(
            set_runtime_activity(Forward), compute_f!, Const,
            Duplicated(p.state, p.dstate),
            Const(p.grid), Const(p.soil), Const(p.constants),
        )
        Jf_vc = vec(interior(p.dstate.tendencies.internal_energy))
        for i in 1:N_z
            for off in -1:1
                j = i + off
                if 1 <= j <= N_z && tridiag_color(j) == c
                    nzv[p.jac_pos[(i, j)]] = ifelse(i == j, one(FT), zero(FT)) - p.Δt * Jf_vc[i]
                end
            end
        end
    end
    return JW
end

# CSC tridiagonal prototype for W, plus a (row, col) → nonzero-position map for fast writes.
W_proto = sparse(Tridiagonal(ones(FT, N_z - 1), ones(FT, N_z), ones(FT, N_z - 1)))

function nonzero_positions(W)
    rv = rowvals(W)
    pos = Dict{Tuple{Int, Int}, Int}()
    for j in axes(W, 2), p in nzrange(W, j)
        pos[(rv[p], j)] = p
    end
    return pos
end

p = merge(p, (; jac_pos = nonzero_positions(W_proto)))

# ## Tracer-based automatic sparsity detection (Route B)
#
# Tracer-based detection through the Float64 state wrapper fails at the `u → state` write
# (`TypeError` from `convert(Float64, tracer)`). What works instead — without touching any
# physics — is re-typing the *container*: build a copy of the `StateVariables` whose
# `Field`s carry SparseConnectivityTracer `Dual` numbers (primal + pattern), then run the
# unmodified `compute_f!` on it. We use the local detector, whose `Dual` tracers carry
# primal values; the global `GradientTracer` has no primal and fails in `safediv`
# (`iszero` requires a primal). Note that comparisons on `Dual` tracers return `Bool`
# (primal-valued), so `ifelse` selects a single branch — the detected pattern reflects the
# branch configuration at the detection state u0; it is verified against the exact
# tridiagonal pattern below.

# Minimal compatibility shim: `SoilComposition`'s inner constructor requires a uniform
# number type for porosity/saturation/liquid/solid, but with a tracer-typed state the
# porosity/solid parameters are Float64 while saturation/liquid are tracers. This
# promoting fallback is scoped to SCT `Dual`s so it never shadows the uniform Float64
# path. Upstream fix candidate: make the `SoilComposition` (and related) constructors
# promote their arguments — required for any non-Enzyme AD through the closures.
function Terrarium.SoilComposition(
        por::Float64, sat::DT, liq::DT, solid::Terrarium.MineralOrganic{Float64}
    ) where {DT <: SparseConnectivityTracer.Dual}
    texture = Terrarium.SoilTexture(DT(solid.texture.sand), DT(solid.texture.clay), DT(solid.texture.silt))
    return Terrarium.SoilComposition(DT(por), sat, liq, Terrarium.MineralOrganic(texture, DT(solid.organic)))
end

# Re-typing helpers: rebuild each `Field` with element type `T`, copying current values
# (including the frozen ghost cells) as constant tracers.
function retype_field(::Type{T}, f::Field{LX, LY, LZ}) where {T, LX, LY, LZ}
    f_new = Field{LX, LY, LZ}(f.grid, T; boundary_conditions = f.boundary_conditions, indices = f.indices)
    f_new.data .= T.(f.data)
    return f_new
end
retype_state(::Type{T}, nt::NamedTuple) where {T} = map(x -> retype_state(T, x), nt)
retype_state(::Type{T}, f::Field) where {T} = retype_field(T, f)
retype_state(::Type{T}, x) where {T} = x
function retype_state(::Type{T}, s::StateVariables) where {T}
    return StateVariables(
        eltype(s), Terrarium.closure_names(s),
        retype_state(T, s.prognostic), retype_state(T, s.tendencies),
        retype_state(T, s.auxiliary), retype_state(T, s.inputs),
        retype_state(T, s.namespaces), s.timestepper_cache, s.clock,
    )
end

sct_detector = TracerLocalSparsityDetector()
TT = jacobian_eltype(zeros(FT, 1), sct_detector)
state_tr = retype_state(TT, state)

# The traced residual: identical formula to `G_nls!`, evaluated on the tracer-typed state.
function G_tr!(res, u)
    u2state!(state_tr, u)
    compute_f!(state_tr, model.grid, model.soil, model.constants)
    res .= u .- u_prev .- Δt .* vec(interior(state_tr.tendencies.internal_energy))
    return nothing
end

@time pattern_sct = jacobian_sparsity(G_tr!, zero(u0), u0, sct_detector)
tridiag_bool = [abs(i - j) <= 1 for i in 1:N_z, j in 1:N_z]
println("SCT-detected pattern: nnz = $(nnz(pattern_sct)) (3N−2 = $(3 * N_z - 2)), " *
        "== tridiagonal: $(collect(pattern_sct .!= 0) == tridiag_bool)")

# ## Dense Enzyme Jacobian (manual Route A detection)
#
# `DenseSparsityDetector(AutoEnzyme)` (variant 4a) cannot differentiate through the
# state-aliasing wrapper (see header note). The same detection semantics are available
# manually: a dense Jacobian from N_z one-hot Enzyme JVPs (the reference example's
# pattern), thresholded to a sparsity pattern. This is the working "Route A with Enzyme".

function dense_jacobian_enzyme(u_eval)
    u2state!(state, u_eval)
    jac = zeros(FT, N_z, N_z)
    shadow_u = interior(p.dstate.prognostic.internal_energy)
    for j in 1:N_z
        shadow_u .= 0
        shadow_u[1, 1, j] = one(FT)
        Enzyme.autodiff(
            set_runtime_activity(Forward), compute_f!, Const,
            Duplicated(state, p.dstate),
            Const(model.grid), Const(model.soil), Const(model.constants),
        )
        jac[:, j] .= vec(interior(p.dstate.tendencies.internal_energy))
    end
    return I(N_z) .- Δt .* jac # W = I − Δt·J_f
end

W_dense = dense_jacobian_enzyme(u0)
offband_dense = maximum(abs.(W_dense .- W_dense .* tridiag_bool))
pattern_enzyme = sparse(map(!iszero, W_dense))
println("dense Enzyme Jacobian: max off-band entry of W = $offband_dense (0 ⟺ exactly tridiagonal), " *
        "pattern == SCT pattern: $(pattern_enzyme == (pattern_sct .!= 0))")

# ## Solver variants

res0 = zero(u0)

nlfunc_jvp = NonlinearFunction(G_nls!; resid_prototype = res0, jvp = jvp_nls!)
nlfunc_jac = NonlinearFunction(G_nls!; resid_prototype = res0, jac = jac_nls!, jac_prototype = W_proto)
nlfunc_autoenzyme = NonlinearFunction(G_nls!; resid_prototype = res0)
nlfunc_route_a = NonlinearFunction(
    G_nls!; resid_prototype = res0,
    sparsity = DenseSparsityDetector(AutoEnzyme(mode = Forward); atol = 0.0),
)
# Route A with Enzyme is blocked by the DI ↔ Enzyme in-place incompatibility (see results);
# finite differences are the working backend for fully automatic sparsity exploitation.
nlfunc_route_a_fd = NonlinearFunction(
    G_nls!; resid_prototype = res0,
    sparsity = DenseSparsityDetector(AutoFiniteDiff(); atol = 0.0),
)
nlfunc_findiff = NonlinearFunction(G_nls!; resid_prototype = res0)
nlfunc_default = NonlinearFunction(G_nls!; resid_prototype = res0)
# Route B: pattern detected by tracers through the unmodified `compute_f!` (above);
# the solve itself runs on the Float64 state with a colored finite-difference Jacobian.
nlfunc_sct = NonlinearFunction(G_nls!; resid_prototype = res0, jac_prototype = pattern_sct)
# Route A with Enzyme, done manually: pattern from the dense Enzyme Jacobian (above).
nlfunc_route_a_manual = NonlinearFunction(G_nls!; resid_prototype = res0, jac_prototype = pattern_enzyme)

variants = [
    (name = "1a. manual Enzyme jvp + GMRES",
        prob = NonlinearProblem(nlfunc_jvp, copy(u0), p),
        alg = NewtonRaphson(linsolve = KrylovJL_GMRES())),
    (name = "1b. manual Enzyme jvp + GMRES + EisenstatWalkerForcing2",
        prob = NonlinearProblem(nlfunc_jvp, copy(u0), p),
        alg = NewtonRaphson(linsolve = KrylovJL_GMRES(), forcing = EisenstatWalkerForcing2())),
    (name = "2. jvp_autodiff=AutoEnzyme + GMRES (experimental)",
        prob = NonlinearProblem(nlfunc_autoenzyme, copy(u0), p),
        alg = NewtonRaphson(linsolve = KrylovJL_GMRES(), jvp_autodiff = AutoEnzyme(mode = Forward))),
    (name = "3. manual colored Enzyme jac + tridiagonal CSC + KLU",
        prob = NonlinearProblem(nlfunc_jac, copy(u0), p),
        alg = NewtonRaphson(linsolve = KLUFactorization())),
    (name = "4a. Route A: DenseSparsityDetector(AutoEnzyme) + KLU (experimental)",
        prob = NonlinearProblem(nlfunc_route_a, copy(u0), p),
        alg = NewtonRaphson(linsolve = KLUFactorization())),
    (name = "4b. automatic sparsity via DenseSparsityDetector(AutoFiniteDiff) + KLU",
        prob = NonlinearProblem(nlfunc_route_a_fd, copy(u0), p),
        alg = NewtonRaphson(linsolve = KLUFactorization(), autodiff = AutoFiniteDiff())),
    (name = "4c. Route A-manual: dense Enzyme jac → pattern + colored FD + KLU",
        prob = NonlinearProblem(nlfunc_route_a_manual, copy(u0), p),
        alg = NewtonRaphson(linsolve = KLUFactorization(), autodiff = AutoFiniteDiff())),
    (name = "5. Route B: TracerLocalSparsityDetector pattern + colored FD + KLU",
        prob = NonlinearProblem(nlfunc_sct, copy(u0), p),
        alg = NewtonRaphson(linsolve = KLUFactorization(), autodiff = AutoFiniteDiff())),
    (name = "6. AutoFiniteDiff dense + LU (baseline)",
        prob = NonlinearProblem(nlfunc_findiff, copy(u0), p),
        alg = NewtonRaphson(autodiff = AutoFiniteDiff())),
    (name = "7. default polyalgorithm (fallback recovery)",
        prob = NonlinearProblem(nlfunc_default, copy(u0), p),
        alg = nothing),
]

results = NamedTuple[]
for v in variants
    local sol
    try
        sol = v.alg === nothing ?
            solve(v.prob; abstol, reltol) :
            solve(v.prob, v.alg; abstol, reltol)
        success = sol.retcode == ReturnCode.Success
        u2state!(state, sol.u)
        compute_f!(state, model.grid, model.soil, model.constants)
        resid = maximum(abs.(sol.u .- u_prev .- Δt .* vec(interior(state.tendencies.internal_energy))))
        println("$(v.name)\n    retcode = $(sol.retcode), |G(u*)|∞ = $resid, stats = $(sol.stats)")
        push!(results, (name = v.name, success = true, u = sol.u, resid = resid, sol = sol))
    catch err
        msg = first(split(sprint(showerror, err), '\n'))
        println("$(v.name)\n    FAILED: $(typeof(err))\n    $msg")
        push!(results, (name = v.name, success = false, u = nothing, resid = nothing, sol = nothing))
    end
end

# ## Ariadne.jl reference solve (Newton–Krylov with finite-difference JVPs)

function f_ariadne!(du, u, p, t)
    state, grid, soil, constants, compute_f! = p
    u2state!(state, u)
    compute_f!(state, grid, soil, constants)
    du .= vec(interior(state.tendencies.internal_energy))
    return nothing
end

function G_euler!(res, u_prev, Δt, f!, du, u, p, t)
    f!(du, u, p, t)
    res .= u .- u_prev .- Δt .* du
    return nothing
end

F_ariadne!(res, u, params) = G_euler!(res, params.u_prev, params.Δt, f_ariadne!, params.du, u, params.p, params.t)

u_ariadne = copy(u_prev)
du_ariadne = zero(u_prev)
res_ariadne = zero(u_prev)
p_ariadne = (state = state, grid = model.grid, soil = model.soil, constants = model.constants, compute_f! = compute_f!)
params_ariadne = (u_prev = u_prev, Δt = Δt, du = du_ariadne, p = p_ariadne, t = FT(0))

_, stats_ariadne = newton_krylov!(F_ariadne!, u_ariadne, params_ariadne, res_ariadne; algo = :gmres, tol_abs = abstol, verbose = 0)
println("Ariadne reference: solved = $(stats_ariadne.solved)")

# ## Correctness checks

# (a) all successful solutions agree with each other and with Ariadne
#     (u is O(10⁷) J/m³, so absolute differences of O(10⁻²) correspond to ~10⁻⁹ relative)
r_success = [r for r in results if r.success]
shortname(name) = first(split(name, '.'))
for r in r_success
    println("agreement |u($(shortname(r.name))) − u(Ariadne)|∞ = $(maximum(abs.(r.u .- u_ariadne)))")
end
for i in eachindex(r_success), j in (i + 1):length(r_success)
    d = maximum(abs.(r_success[i].u .- r_success[j].u))
    println("agreement |u($(shortname(r_success[i].name))) − u($(shortname(r_success[j].name)))|∞ = $d")
end

# (b) the colored decoding of `jac_nls!` is checked against the dense Enzyme reference
#     Jacobian computed above (both at u0); tridiagonality was verified there
W_colored = copy(W_proto)
jac_nls!(W_colored, u0, p)
println("colored Enzyme jac vs dense Enzyme reference: max |ΔW| = $(maximum(abs.(W_colored - W_dense)))")

# compact summary
println("\n=== summary ===")
for r in results
    status = r.success ? "|G(u*)|∞ = $(r.resid)" : "FAILED"
    println("$(rpad(shortname(r.name), 4)) $(status)")
end

# ## Benchmarks
#
# Each successful variant is benchmarked on the same problem. The first solve above
# triggered compilation, so these timings reflect steady-state performance.

for (i, v) in enumerate(variants)
    results[i].success || continue
    println("benchmark: $(v.name)")
    prob_v, alg_v = v.prob, v.alg
    if alg_v === nothing
        display(@benchmark solve($prob_v; abstol = $abstol, reltol = $reltol))
    else
        display(@benchmark solve($prob_v, $alg_v; abstol = $abstol, reltol = $reltol))
    end
end

# ## Stage 2 design notes (kernel-launched per-column solves)
#
# Stage 2 lifts this single-column experiment to a full grid: one N_z-dimensional
# nonlinear solve per soil column, all columns solved concurrently by launching the
# Newton iteration itself as a KernelAbstractions kernel over columns. Columns are
# independent (1D vertical physics), so the global Jacobian is block-diagonal and no
# cross-column communication is needed.
#
# Consequences for the algorithm choice:
#   - Krylov methods (GMRES) are not viable inside a kernel: no ecosystem support, and
#     the O(N_z) tridiagonal structure makes them unnecessary anyway.
#   - The candidate in-kernel solvers are matrix-free quasi-Newton methods from
#     SimpleNonlinearSolve (`SimpleKlement`, `SimpleLimitedMemoryBroyden`, `SimpleDFSane`)
#     or a custom tridiagonal Newton with a Thomas-algorithm solve. StaticArrays `\` is
#     unreliable for N ≳ 14, so the Thomas solve should be hand-written.
#   - The residual inside the kernel must be allocation-free, throw-free (no reachable
#     trap paths), and use `ifelse` instead of branches — per the Terrarium kernel rules.
#     The tracer experiment above shows that the existing kernel functions
#     (`energy_to_temperature!`, `compute_energy_tendency`, …) are themselves
#     eltype-generic, so the in-kernel residual can be assembled from Terrarium's own
#     kernel functions (called per-index over plain buffers) rather than re-coded; the
#     main blocker found is the uniform-number-type constructor gate in `SoilComposition`
#     (see the shim above).
#   - The Newton matrix W = I − Δt·J_f is tridiagonal; its bands have closed-form
#     expressions through the (piecewise-linear in liq) freeze curve and the
#     (A + B·liq)² conductivity form, so an analytic in-kernel Jacobian is feasible
#     without any AD.
