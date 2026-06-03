# # Differentiating Terrarium.jl
#
# We build Terrarium with differentiability in mind. This means that you are able to take derivatives of outputs of Terrarium with automatic differentiation (AD). AD enables us to use e.g. an automated, objetive calibration of model parameters, but also the direct integration of neural networks and other machine learning (ML) methods into our model. For this purpose we ensure compatibility on [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl). Enzyme.jl can peform both a reverse-mode AD (typical for most ML applications), and a forward-mode AD (more typical in classical sensitivity analysis).
#
# When differentiating through a model integration, AD would usually need to keep track of every single intermediate value that contributes to our output. For long integrations this quickly becomes infeasible due to its high memory demand. Therefore we support checkpointing schemes from [Checkpointing.jl](https://github.com/Argonne-National-Laboratory/Checkpointing.jl) for these cases that only save selected intermediate values and recompute all other intermediate values when needed. For Enzyme.jl this also has another very practical advantage: the first compile time for the gradient is much lower.
#
# Without much further ado, let us look into how we can differentiate Terrarium hands-on and perform a small sensitivity analysis of a one column soil model. First, we set up our model as usual:

using Terrarium, Enzyme, Checkpointing
using LinearAlgebra

import CairoMakie as Makie

grid = ColumnGrid(UniformSpacing()) # Easier Jacobian
initializer = SoilInitializer(eltype(grid))
model = SoilModel(grid; initializer)
# constant surface temperature of 1°C
bcs = PrescribedSurfaceTemperature(:T_ub, 1.0)
integrator = initialize(model, ForwardEuler(), boundary_conditions = bcs)

# So far, this is just our usual setup. In this case, for a soil column with a prescribed surface temperature.
#
# Now, we set up our AD checkpoiting scheme for the timestepping. Here we choose a [Revolve scheme](https://dl.acm.org/doi/10.1145/347837.347846) that saves intermediate values at every single time step. Note that when we save at every single time step the different available schemes don't actually differ from each other.

scheme = Revolve(1)

# Next we prepare to differentiate with Enzyme. For a comprehensive introduction to Enzyme, please see [their documentation](https://enzymead.github.io/Enzyme.jl/stable/).
#
# We want to perform a sensitivity analysis of the temperature of the second lowest soiler layer ``T_f`` at the end of our simulation with respect to the initial conditions of our simulation ``\mathbf{U}_0``, ``\mathbf{T}_0``, where ``\mathbf{U}`` is the internal energy.
#
# Enzyme's `autodiff` is it's core function that we can use to compute vector-Jacobian products (vJP) for the reverse-mode AD of our `run!` function that integrates our model using the `integrator` that we initialized. In order to compute the gradient of the just one layer of the soil, we set a "one-hot" seed for the vJP like so:

dintegrator = make_zero(integrator)
# set a one hot seed for a sensitivity analysis of T for now
interior(dintegrator.state.temperature)[1, 1, 2] = 1.0

# how many steps we want to integrate for
N_t = 200

# While doing that we allocated a shadow memory `dintegrator` for Enzyme in which it can accumluate the vJP (see Enzyme docs for more information). That's all the setup we need, as we have a pre-defined version of `run!` that takes in our `scheme`. We just need to call `autodiff` now. Executing this for the first time, might take a few minutes. Subsequent executions will be very fast though.

autodiff(set_runtime_activity(Reverse), run!, Const, Duplicated(integrator, dintegrator), Const(scheme), Const(N_t))

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

# As expected the sensitivity is the highest locally, with the same and neighbouring soil layers contributing and no sensitivity wrt higher soil layers for our still rather short integration of only ``N_t\cdot 300s``.
#
# This example should just demonstrate the technical possibilities of Terrarium.jl in an easy and fast to compute setup, stay tuned for more complex examples.

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

N_z = size(interior(state.internal_energy), 3)

# Combined function: closure + tendency in one call so the full U→T→∂U∂t chain
# is visible to Enzyme.
function compute_f!(state, grid, soil, constants)
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

for (col_offset, var_in) in enumerate(prog_vars)
    for k in 1:N_z
        col = (col_offset - 1) * N_z + k

        dstate = make_zero(state)
        interior(getproperty(dstate.prognostic, var_in))[1, 1, k] = one(eltype(grid))

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
            jac_full[rows, col] .= interior(getproperty(dstate.tendencies, var_out))[1, 1, :]
        end
    end
end

Δt = integrator.timestepper.Δt   # or any chosen implicit timestep
W = I(N_z) - Δt * jac_full

f3 = Makie.Figure(size = (800, 400))
ax1 = Makie.Axis(f3[1, 1], yreversed = true, title = "J_f = ∂(∂U∂t)/∂U  (heat equation, NoFlow)", xlabel = "Layer k", ylabel = "Layer i")
ax2 = Makie.Axis(f3[1, 2], yreversed = true, title = "W = I − Δt·J_f  (Newton matrix, Hairer & Wanner)", xlabel = "Layer k", ylabel = "Layer i")
Makie.heatmap!(ax1, jac_full)
Makie.heatmap!(ax2, W)
f3
