# # Differentiating Terrarium.jl
#
# We build Terrarium with differentiability in mind. This means that you are able to take derivatives of outputs of Terrarium with automatic differentiation (AD). AD enables us to use e.g. an automated, objetive calibration of model parameters, but also the direct integration of neural networks and other machine learning (ML) methods into our model. For this purpose we ensure compatibility on [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl). Enzyme.jl can peform both a reverse-mode AD (typical for most ML applications), and a forward-mode AD (more typical in classical sensitivity analysis).
#
# When differentiating through a model integration, AD would usually need to keep track of every single intermediate value that contributes to our output. For long integrations this quickly becomes infeasible due to its high memory demand. Therefore we support checkpointing schemes from [Checkpointing.jl](https://github.com/Argonne-National-Laboratory/Checkpointing.jl) for these cases that only save selected intermediate values and recompute all other intermediate values when needed. For Enzyme.jl this also has another very practical advantage: the first compile time for the gradient is much lower.
#
# Without much further ado, let us look into how we can differentiate Terrarium hands-on and perform a small sensitivity analysis of a one column soil model. First, we set up our model as usual:

using Terrarium
using Enzyme
using Enzyme: Forward, Reverse, set_runtime_activity
using Checkpointing
using CUDA
using BenchmarkTools
using LinearAlgebra
using CairoMakie
using Statistics

FT = Float64
grid = ColumnGrid(CPU(), FT, UniformSpacing(), 1) # Easier Jacobian
initializer = SoilInitializer(eltype(grid))
model = SoilModel(grid; initializer)
# constant surface temperature of 1°C
bcs = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator = initialize(model, boundary_conditions = bcs)

# So far, this is just our usual setup. In this case, for a soil column with a prescribed surface temperature.
#
# Now, we set up our AD checkpointing scheme for the timestepping. Here we choose a [Revolve scheme](https://dl.acm.org/doi/10.1145/347837.347846) that saves intermediate values at every single time step. Note that when we save at every single time step the different available schemes don't actually differ from each other.

scheme = Revolve(1)

# Next we prepare to differentiate with Enzyme. For a comprehensive introduction to Enzyme, please see [their documentation](https://enzymead.github.io/Enzyme.jl/stable/).
#
# We want to perform a sensitivity analysis of the temperature of the second lowest soil layer ``T_f`` at the end of our simulation with respect to the initial conditions of our simulation ``\mathbf{U}_0``, ``\mathbf{T}_0``, where ``\mathbf{U}`` is the internal energy.
#
# Enzyme's `autodiff` is it's core function that we can use to compute vector-Jacobian products (vJP) for the reverse-mode AD of our `run!` function that integrates our model using the `integrator` that we initialized. In order to compute the gradient of the just one layer of the soil, we set a "one-hot" seed for the vJP like so:

dintegrator = make_zero(integrator)
# set a one hot seed for a sensitivity analysis of T for now
interior(dintegrator.state.temperature)[1, 1, 2] = FT(1.0)

function mean_temperature(integrator)
    run!(integrator, steps = 1)
    return mean(interior(integrator.state.temperature))
end

# While doing that we allocated a shadow memory `dintegrator` for Enzyme in which it can accumulate the vJP (see Enzyme docs for more information).
# We just need to call `autodiff` now. Executing this for the first time, might take a few minutes. Subsequent executions will be very fast though.

autodiff(set_runtime_activity(Reverse), mean_temperature, Active, Duplicated(integrator, dintegrator))

# Let's look at the results that were accumulated in our shadow memory `dintegrator` by Enzyme and plot them!

dU = interior(dintegrator.state.internal_energy)[1, 1, :]
dT = interior(dintegrator.state.temperature)[1, 1, :]
zs = znodes(integrator.state.temperature)

f = Makie.Figure()
Makie.Axis(f[1, 1], ylabel = "Soil depth", xlabel = "Sensitivity dT_f/dU_0")
Makie.scatterlines!(f[1, 1], dT, zs)
f

f2 = Makie.Figure()
Makie.Axis(f2[1, 1], ylabel = "Soil depth", xlabel = "Sensitivity dU_f/dU_0")
Makie.scatterlines!(f2[1, 1], dU, zs)
f2

# As expected the sensitivity is the highest locally, with the same and neighboring soil layers contributing and no sensitivity wrt higher soil layers for our still rather short integration of only ``N_t\cdot 300s``.

# ## Parameter sensitivities
#
# We can also compute sensitivities with respect to *model parameters* rather than initial
# conditions. Since we are interested in only a single parameter — the soil mineral thermal
# conductivity ``\lambda_\text{mineral}`` — forward-mode AD is the natural choice: it
# propagates one tangent vector forward through the computation in a single pass, without
# the memory overhead of reverse mode.
#
# We re-initialize a fresh integrator so the state starts at ``t = 0``:

grid = ColumnGrid(ExponentialSpacing())
initializer = SoilInitializer(eltype(grid))
model = SoilModel(grid; initializer)
bcs = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator = initialize(model, ForwardEuler(), boundary_conditions = bcs)

# TODO

function mean_temperature(clock, model, inputs, state, inits, timestepper)
    integrator = Terrarium.ModelIntegrator(clock, model, inputs, state, inits, timestepper)
    run!(integrator, steps = 1)
    return mean(interior(integrator.state.temperature))
end

scheme = Revolve(1)
dmodel = Ref(make_zero(model))
dintegrator = make_zero(integrator)
dstate = dintegrator.state

grads = autodiff(
    set_runtime_activity(Reverse), mean_temperature,
    Active,
    Const(integrator.clock),
    MixedDuplicated(model, dmodel),
    Const(integrator.inputs),
    Duplicated(integrator.state, dstate),
    Const(integrator.initializers),
    Duplicated(integrator.timestepper, make_zero(integrator.timestepper))
)

dTdp = dstate.temperature ./ dmodel.x.soil.energy.thermal_properties.conductivities.mineral
zs = znodes(integrator.state.temperature)
lines(dTdp[1, 1, :])

# The temperature tangent accumulated in `dintegrator` now holds
# ``\partial T(k) / \partial \lambda_\text{mineral}`` for each layer ``k``:

dT_mineral = interior(dintegrator.state.temperature)[1, 1, :]
zs = znodes(integrator.state.temperature)

f3 = Makie.Figure()
Makie.Axis(f3[1, 1], ylabel = "Soil depth", xlabel = "Sensitivity ∂T/∂λ_mineral  (K m K W⁻¹)")
Makie.scatterlines!(f3[1, 1], dT_mineral, zs)
f3

# A positive sensitivity in the upper layers is consistent with a more conductive mineral
# matrix transferring surface warmth downward more efficiently over one time step.
#
# This example should just demonstrate the technical possibilities of Terrarium.jl in an
# easy and fast to compute setup, stay tuned for more complex examples.

# ## Jacobian of the tendency map w.r.t. all prognostic variables
#
# The motivation is **implicit timestepping** of the coupled land model.  A backward-Euler
# step requires solving the nonlinear system
#
#   g(U^{n+1}) = U^{n+1} - U^n - Δt·f(U^{n+1}) = 0
#
# via Newton iterations, each of which needs the Jacobian
#
#   W = I - Δt·J_f,   J_f = ∂f/∂U   (W is the standard notation in Hairer & Wanner)
#
# where U is the prognostic vector: [internal_energy] only in this example.
# The model uses the default `SoilHydrology` with `NoFlow` vertical flow (no Richards
# equation), so `saturation_water_ice` is frozen in place — its tendency is purely
# from freeze/thaw driven by the energy balance and is zero here.
# We therefore build the Jacobian only w.r.t. `internal_energy`.
#
# Key insight: the tendency f(U) flows through two stages
#
#   U  ──closure!──►  T, liq, ...  ──compute_tendencies!──►  ∂U∂t
#
# `compute_auxiliary!(state, grid, energy::SoilEnergyBalance, ...)` is a no-op; the
# actual U→T mapping lives in `closure!` (called by the timestepper after each step).
# Differentiating only `compute_tendencies!` gives zero because temperature hasn't been
# updated from the perturbed energy.  We must differentiate through `closure!` + tendencies.

Terrarium.initialize!(integrator)
state = integrator.state
Δt = integrator.timestepper.Δt   # or any chosen implicit timestep

N_z = size(interior(state.internal_energy), 3)

# Combined function: closure + tendency in one call so the full U→T→∂U∂t chain
# is visible to Enzyme.
function compute_f!(state, grid, soil, constants)
    Terrarium.reset_tendencies!(state)
    Terrarium.closure!(state, grid, soil, constants)
    Terrarium.compute_tendencies!(state, grid, soil, constants)
    return nothing
end

# Only internal_energy has a non-trivial tendency with NoFlow hydrology.
# saturation_water_ice is immobile (NoFlow); its tendency is zero in this configuration.
prog_vars = (:internal_energy,)
n_prog = length(prog_vars) * N_z        # total prognostic DOFs
jac_full = zeros(eltype(grid), n_prog, n_prog)

tend_vars = (:internal_energy,)

# Benchmark a JVP
dstate = make_zero(state)
CUDA.@allowscalar  interior(getproperty(dstate.prognostic, :internal_energy))[1, 1, 1] = one(eltype(grid))

# Measure time for copmpilation
@time Enzyme.autodiff(
    set_runtime_activity(Forward),
    compute_f!,
    Const,
    Duplicated(state, dstate),
    Const(model.grid),
    Const(model.soil),
    Const(model.constants),
)

# Benchmark after compilation for the first time
dstate = make_zero(state)
CUDA.@allowscalar  interior(getproperty(dstate.prognostic, :internal_energy))[1, 1, 1] = one(eltype(grid))

@benchmark Enzyme.autodiff(
    set_runtime_activity(Forward),
    compute_f!,
    Const,
    Duplicated(state, dstate),
    Const(model.grid),
    Const(model.soil),
    Const(model.constants),
)

# Jacobian
for (col_offset, var_in) in enumerate(prog_vars)
    for k in 1:N_z
        col = (col_offset - 1) * N_z + k

        dstate = make_zero(state)
        CUDA.@allowscalar  interior(getproperty(dstate.prognostic, var_in))[1, 1, k] = one(eltype(grid))

        Enzyme.autodiff(
            set_runtime_activity(Forward),
            compute_f!,
            Const,
            Duplicated(state, dstate),
            Const(model.grid),
            Const(model.soil),
            Const(model.constants),
        )

        for (row_offset, var_out) in enumerate(tend_vars)
            rows = ((row_offset - 1) * N_z + 1):(row_offset * N_z)
            CUDA.@allowscalar jac_full[rows, col] .= interior(getproperty(dstate.tendencies, var_out))[1, 1, :]
        end
    end
end

W = I(N_z) - Δt * jac_full

f3 = Makie.Figure(size = (800, 400))
ax1 = Makie.Axis(f3[1, 1], yreversed = true, title = "J_f = ∂(∂U∂t)/∂U  (heat equation, NoFlow)", xlabel = "Layer k", ylabel = "Layer i")
ax2 = Makie.Axis(f3[1, 2], yreversed = true, title = "W = I − Δt·J_f  (Newton matrix, Hairer & Wanner)", xlabel = "Layer k", ylabel = "Layer i")
Makie.heatmap!(ax1, jac_full)
Makie.heatmap!(ax2, W)
f3

# Using Newton-Krylov method
using Ariadne

# Basically just copying from https://numericalmathematics.github.io/Ariadne.jl/stable/generated/implicit/
# Renaming u_n to u_prev for clarity: u_prev is the known state at the previous time step.

# 1. Adapt your compute_f! to the Ariadne interface: f!(du, u, p, t)
# This wrapper handles the mapping between the vector `u` used by the solver
# and the structured `state` used by Terrarium.
function f_ariadne!(du, u, p, t)
    state, grid, soil, constants, compute_f! = p

    # Update internal state with current guess u
    interior(state.prognostic.internal_energy) .= reshape(u, size(interior(state.prognostic.internal_energy))) # back to 3D shape

    # Compute the tendencies f(u)
    compute_f!(state, grid, soil, constants)

    # Extract the computed tendency into du
    du .= vec(interior(state.tendencies.internal_energy)) # 3D to 1D shape
    return nothing
end

# 2. Define the Implicit Euler residual: G = u - u_prev - Δt*f(u) = 0
function G_euler!(res, u_prev, Δt, f!, du, u, p, t)
    f!(du, u, p, t)
    res .= u .- u_prev .- Δt .* du
    return nothing
end

function jacobian(G!, f!, u_prev, p, Δt, t)
    u = copy(u_prev)
    du = zero(u_prev)
    res = zero(u_prev)

    F!(res, u, (u_prev, Δt, du, p, t)) = G!(res, u_prev, Δt, f!, du, u, p, t)

    J = Ariadne.JacobianOperator(F!, res, u, (u_prev, Δt, du, p, t))
    return collect(J)
end

# 3. Setup the solve

Terrarium.initialize!(integrator) # reset integrator
state = integrator.state
u_prev = copy(vec(interior(state.prognostic.internal_energy)))
u = copy(u_prev)      # Initial guess (usually the solution at the previous step)
du = zero(u)
res = zero(u)
t = FT(0.0)        # Time (not used in this autonomous example, but required for signature)

# Parameters passed into f_ariadne!
p = (state = state, grid = model.grid, soil = model.soil, constants = model.constants, compute_f! = compute_f!)

# Define the wrapper F!(res, u, params) that Ariadne's solver expects
# The params tuple matches the arguments of G_euler!
F!(res, u, params) = G_euler!(res, params.u_prev, params.Δt, f_ariadne!, params.du, u, params.p, params.t)

params = (u_prev = u_prev, Δt = Δt, du = du, p = p, t = t)

# 4. Perform the Newton-Krylov solve
# Note: For energy units (10^7 J/m³), we use a more realistic absolute tolerance.
newton_krylov_kwargs = (tol_abs = 1.0e-2, verbose = 0)
algo = :gmres

@time  _, stats = newton_krylov!(F!, u, params, res; algo, newton_krylov_kwargs...)

if stats.solved
    println("Newton-Krylov converged")
    # The converged u now contains the updated internal_energy
else
    @error "Newton-Krylov solve failed!" stats
end

# 5. Optional: Compute the Jacobian accurately using Ariadne's JacobianOperator
# This can be used to compare against your manual jac_full
J_op = Ariadne.JacobianOperator(F!, res, u, params)
jac_ariadne = collect(J_op)

fig = Makie.Figure(size = (400, 400))
ax = Makie.Axis(fig[1, 1], yreversed = true, title = "W Jacobian from Ariadne")
Makie.heatmap!(ax, jac_ariadne)
fig
