# # Differentiating Terrarium.jl
#
# We build Terrarium with differentiability in mind. This means that you are able to take derivatives of outputs of Terrarium with automatic differentiation (AD). AD enables us to use e.g. an automated, objetive calibration of model parameters, but also the direct integration of neural networks and other machine learning (ML) methods into our model. For this purpose we ensure compatibility on [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl). Enzyme.jl can peform both a reverse-mode AD (typical for most ML applications), and a forward-mode AD (more typical in classical sensitivity analysis).
#
# When differentiating through a model integration, AD would usually need to keep track of every single intermediate value that contributes to our output. For long integrations this quickly becomes infeasible due to its high memory demand. Therefore we support checkpointing schemes from [Checkpointing.jl](https://github.com/Argonne-National-Laboratory/Checkpointing.jl) for these cases that only save selected intermediate values and recompute all other intermediate values when needed. For Enzyme.jl this also has another very practical advantage: the first compile time for the gradient is much lower.
#
# Without much further ado, let us look into how we can differentiate Terrarium hands-on and perform a small sensitivity analysis of a one column soil model. First, we set up our model as usual:

using Terrarium, Enzyme, Checkpointing

import CairoMakie as Makie

grid = ColumnGrid(ExponentialSpacing())
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
