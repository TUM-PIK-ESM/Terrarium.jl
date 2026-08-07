# # Differentiating Terrarium.jl
#
# We build Terrarium with differentiability in mind. This means that you are able to take derivatives of outputs of Terrarium with automatic differentiation (AD). AD enables us to use e.g. an automated, objective calibration of model parameters, but also the direct integration of neural networks and other machine learning (ML) methods into our model. For this purpose we ensure compatibility on [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl). Enzyme.jl can perform both a reverse-mode AD (typical for most ML applications), and a forward-mode AD (more typical in classical sensitivity analysis).
#
# When differentiating through a model integration, AD would usually need to keep track of every single intermediate value that contributes to our output. For long integrations this quickly becomes infeasible due to its high memory demand. Therefore we support checkpointing schemes from [Checkpointing.jl](https://github.com/Argonne-National-Laboratory/Checkpointing.jl) for these cases that only save selected intermediate values and recompute all other intermediate values when needed. For Enzyme.jl this also has another very practical advantage: the first compile time for the gradient is much lower.
#
# Without much further ado, let us look into how we can differentiate Terrarium hands-on and perform a small sensitivity analysis of a one column soil model. First, we set up our model as usual:

using Terrarium
using Enzyme
using Enzyme: Forward, Reverse, set_runtime_activity
using Checkpointing
using CairoMakie
using Statistics

arch = CPU()
FT = Float32
grid = ColumnGrid(arch, FT, ExponentialSpacing())
initializer = SoilInitializer(FT)
model = SoilModel(grid; timestepper = ForwardEuler(FT), initializer = initializer)
# constant surface temperature of 1°C
bcs = PrescribedSurfaceTemperature(:T_ub, FT(1.0))
integrator = initialize(model, boundary_conditions = bcs)

# So far, this is just our usual setup. In this case, for a soil column with a prescribed surface temperature.
#
# Now, we set up our AD checkpointing scheme for the timestepping. Here we choose a [Revolve scheme](https://dl.acm.org/doi/10.1145/347837.347846) that saves intermediate values at every single time step. Note that when we save at every single time step the different available schemes don't actually differ from each other.

scheme = Revolve(1)

# Next we prepare to differentiate with Enzyme. For a comprehensive introduction to Enzyme, please see [their documentation](https://enzymead.github.io/Enzyme.jl/stable/).
#
# We want to perform a sensitivity analysis of the temperature of the second lowest soil layer ``T_f`` at the end of our simulation with respect to the initial conditions of our simulation ``\mathbf{U}_0``, ``\mathbf{T}_0``, where ``\mathbf{U}`` is the internal energy.
#
# Enzyme's `autodiff` is it's core function that we can use to compute vector-Jacobian products (vJP) for the reverse-mode AD of our `run!` function that integrates our model using the `integrator` that we initialized. We define a function that returns the temperature of the second lowest soil layer after a time integration of 10 steps.

dintegrator = make_zero(integrator)

function layer_temperature(integrator)
    run!(integrator, scheme, 10)
    return interior(integrator.state.temperature)[1, 1, 2]
end

# While doing that we allocated a shadow memory `dintegrator` for Enzyme in which it can accumulate the vJP (see Enzyme docs for more information).
# We just need to call `autodiff` now. Some notes ion its arguments:
# - `set_runtime_activity(Reverse)` tells Enzyme to use reverse-mode AD while enabling runtime activity (see [here](https://enzymead.github.io/Enzyme.jl/dev/faq/#faq-runtime-activity) for details)
# - `Active` tells Enzyme that we want to take the derivative of the scalar output of `layer_temperature` with respect to the function's input. The cotangent of the function output is set to 1.0 in this way for the vJP calculation.
# - `Duplicated(integrator, dintegrator)` `dintegrator` is shadow memory that Enzyme uses to accumulate the vJP of the `integrator` state variables.

#Executing this for the first time, might take a few minutes. Subsequent executions will be very fast though.

autodiff(set_runtime_activity(Reverse), layer_temperature, Active, Duplicated(integrator, dintegrator))

# Let's look at the results that were accumulated in our shadow memory `dintegrator` by Enzyme and plot them!

dU = interior(dintegrator.state.internal_energy)[1, 1, :]
dT = interior(dintegrator.state.temperature)[1, 1, :]
zs = znodes(integrator.state.temperature)

f = Makie.Figure()
Makie.Axis(f[1, 1], ylabel = "Soil depth", xlabel = "Sensitivity dT_f/dT_0")
Makie.scatterlines!(f[1, 1], dT, zs)
f

f2 = Makie.Figure()
Makie.Axis(f2[1, 1], ylabel = "Soil depth", xlabel = "Sensitivity dT_f/dU_0")
Makie.scatterlines!(f2[1, 1], dU, zs)
f2

# As expected the sensitivity is the highest locally, with the same and neighboring soil layers contributing and no sensitivity wrt higher soil layers for our still rather short integration.

# ## Parameter sensitivities
#
# We can also compute sensitivities with respect to *model parameters* rather than initial
# conditions. Since we are interested in only a single parameter — the soil mineral thermal
# conductivity ``\lambda_\text{mineral}`` — forward-mode AD is the natural choice: it
# propagates one tangent vector forward through the computation in a single pass, without
# the memory overhead of reverse mode.
#
# We re-initialize a fresh integrator so the state starts at ``t = 0``:

grid = ColumnGrid(arch, FT, ExponentialSpacing())
initializer = SoilInitializer(FT)
model = SoilModel(grid; timestepper = ForwardEuler(FT), initializer = initializer)
bcs = PrescribedSurfaceTemperature(:T_ub, FT(1.0))
integrator = initialize(model, boundary_conditions = bcs)

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
