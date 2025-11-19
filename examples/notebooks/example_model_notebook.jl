### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 808d5d89-c1d2-4f6a-bd55-4b3a8444c90f
begin
    import Pkg
    Pkg.activate(".")
end

# ╔═╡ 94d82d31-42ec-41de-91e9-b5585b3a72d4
using Terrarium

# ╔═╡ 07c8a3a4-21aa-4213-a876-eadc8754d4a0
begin # for reproducability set a random seed
	using Random
	Random.seed!(1234)
end

# ╔═╡ eeb283fa-5360-4bab-83cf-dcbc0bee7949
using CairoMakie

# ╔═╡ 5630efd5-2482-463d-913f-9addb120beec
md"""
# A super basic example model with Terrarium 

In this example we will set up an embarrassingly simple example to demonstrate Terrarium's model interface. Our model will have 1-dimensional exponential dynamics with a constant offset

```math 
\frac{du}{dt} = \alpha u + c + F(t)
```

for an arbitrary prognostic variable ``u``. For the sake of this demonstration we will treat the offset ``c`` as an auxiliary/diagnostic variable even though it is constant in time. ``F(t)`` is an external forcing that we apply. 

We begin by defining our model `struct` that subtypes `Terrarium.AbstractModel`: 
"""

# ╔═╡ 4922e264-c80d-4a5b-8891-a7c8a3fdbfe7
md""" 
A "model" in Terrarium is a subtype of `Terrarium.AbstractModel` and is a `struct` type constisting of 
 * `grid` which defines the discretization of the spatial domain
 * `initializer` which is responsible for initializing state variables
 * further fields that define processes, dynamics and submodels 

When we follow the advised naming notations of `grid` and `initializer` we inherit default methods from `Terrarium.AbstractModel` such as `get_grid` and `get_initializer`. For more complex models we might need to implement custom overrides of `initialize!(state, ::Model, ::Initializer)` to initialize model states. 

## What is a "grid"? 

The `grid` defines the spatial discretization. Our grids are based on those of Oceananigans.jl (and SpeedyWeather.jl/RingGrids.jl) in order to take advantage of their capabilities for device-agnostic parallelization.

Terrarium currently provides two grid types: 

* `ColumnGrid` is a set of laterally independent vertical columns with dimensions ``(x, y, z)`` where ``x`` is the column dimension, ``y=1`` is constant, and ``z`` is the vertical axis, 
* `ColumnRingGrid` represents a global (spherical) grid of independent, vertical columns where the spatial discretization in the horizontal direction is defined by a `RingGrids.AbstractGrid`. 

In both cases we need to specificy the vertical discretizataion via an `UniformSpacing`, `ExponentialSpacing` or `PrescribedSpacing`.

## Initializer and Boundary Conditions

For our basic example here the default initializer (which does nothing) will suffice, and we won't have to define a custom one.

Boundary conditions are specified by passing Oceananigans `BoundaryCondition` types to `initialize`. In the case of a linear ODE, however, no boundary conditions are required.

## What's our `grid`? 

For our current example, we are defining a simple linear ODE without any spatial dynamics, so we can get away with just a single column with one vertical layer. We can define it like so:
"""

# ╔═╡ 78f268ef-5385-4c63-bc35-2c973de69da5
grid = ColumnGrid(CPU(), Float64, UniformSpacing(N=1))

# ╔═╡ 054a8b11-250f-429f-966f-ca3c9a5dc2ef
md"""
## Defining the model

We start by defining a `struct` for our model that inherits from `AbstractModel` and consists of three properties: the spatial `grid`, an `initializer`, and a single `AbstractProcess` defining the dynamics, which we will also implement below.
"""

# ╔═╡ 407786db-607f-4508-b697-fe75b3ce0b25
begin
    @kwdef struct LinearDynamics{NF} <: Terrarium.AbstractProcess
        "Exponential growth rate"
        alpha::NF = 0.01

        "Constant offset for exponential growth"
        c::NF = 0.1
    end

    @kwdef struct ExpModel{NF, Grid <: Terrarium.AbstractLandGrid{NF}, Dyn, Init} <: Terrarium.AbstractModel{NF, Grid}
        "Spatial grid on which state variables are discretized"
        grid::Grid

        "Linear dynamics resulting in exponential growth/decay"
        dynamics::Dyn = LinearDynamics()

        "Model initializer"
        initializer::Init = DefaultInitializer()
    end
end

# ╔═╡ 575d920c-b12e-493f-95a7-5c962c3591fd
md"""
## Defining the model behaviour 

Now, we want to define our intended model behaviour. For this, we need to define the following methods: 

* `variables(::Model)` returns a tuple of variable metadata declaring the state variables. Variables must be one of three types: `prognostic`, `auxiliary` (sometimes referred to as “diagnostic”), or `input`. Prognostic variables fully characterize the state of the system at any given timestep and are updated according to their tendencies (i.e. ``u`` in our example). Tendencies are automatically allocated for each prognostic variable declared by the model. In this example we will treat the offset ``c`` as an auxiliary variable, though we could also just include it as a constant in the tendency computations.
* `compute_auxiliary!(state, ::Model)` computes the values of all auxiliary variables (if necessary) assuming that the prognostic variables of the system in state are available for the current timestep.
* `compute_tendencies!(state, ::Model)` computes the tendencies based on the current values of the prognostic and auxiliary variables stored in state.

So, let's define those: 
"""

# ╔═╡ 82e45724-ba16-4806-9470-5cb4c43ea734
Terrarium.variables(::ExpModel) = (
    Terrarium.prognostic(:u, Terrarium.XY()),
    Terrarium.auxiliary(:c, Terrarium.XY()),
	Terrarium.input(:F, Terrarium.XY()),
)

# ╔═╡ d4d19de7-6f77-4873-9182-9832d1ca4381
md"""
Here, we defined our three variables with their names as a `Symbol` and whether they are 2D variables (`XY`) on the spatial grid or 3D variables (`XYZ`) that also vary along the vertical z-axis. Here we are considering only a simple scalar model so we choose 2D (`XY`), bearing in mind that all points in the X and Y dimensions of `ColumnGrid` are independent of each other.

We also need to define `compute_auxiliary!` and `compute_tendencies!` as discussed above. We will use here a pattern which is commonly employed within Terrarium: we unpack the process from the model and forward the method calls to more specialzied ones defined for the `LinearDynamics` process.
"""

# ╔═╡ 5ea313fc-3fbb-4092-a2cc-e0cd1f2fe641
function Terrarium.compute_auxiliary!(state, model::ExpModel)
    compute_auxiliary!(state, model, model.dynamics)
end

# ╔═╡ 3815424f-6210-470d-aef1-99c60c71072f
function Terrarium.compute_tendencies!(state, model::ExpModel)
    compute_tendencies!(state, model, model.dynamics)
end

# ╔═╡ 32373599-768f-4809-acdd-4704acc3f30b
md"""
Note that, when implementing models within the Terrarium module itself, the `Terrarium.` qualifier in the definition is not needed.

## Implementing the dynamics

Next, we define the functions that compute the actual dynamics. In order to do this, we need to know a little about how the variables we just defined are handled in our `StateVariables`. The `StateVariables` hold all prognostic and auxiliary variables, their tendencies and closures and additional inputs and forcings in seperate `NamedTuples`. Note that Terrarium also defines shortcuts such that, e.g. in our example, both `state.prognostic.u` and `state.u` would work.

With that in mind, let's define the methods:
"""

# ╔═╡ d55aaf4c-3033-45ba-9d64-8fa8ae4b671c
function Terrarium.compute_auxiliary!(
        state,
        model::ExpModel,
        dynamics::LinearDynamics
    )
    # set auxiliary variable for offset c
    return state.auxiliary.c .= dynamics.c
end

# ╔═╡ 5c8be7e4-f150-492b-a75d-96887a11f6da
# du/dt = u + c
function Terrarium.compute_tendencies!(
        state,
        model::ExpModel,
        dynamics::LinearDynamics
    )
    # define the dynamics; we'll use some special characters to make the equation nicer to look at :)
    return let u = state.prognostic.u,
            ∂u∂t = state.tendencies.u,
            α = dynamics.alpha,
            # note again that here we could also just use dynamics.c instead of defining an auxiliary variable!
            c = state.auxiliary.c
			F = state.inputs.F
        # Write into tendency variable ∂u∂t
        ∂u∂t .= α * u + c + F
    end
end

# ╔═╡ 8d331856-6e9b-41d4-b1be-a84d5fedac8d
md"""
These example compute functions are really the simplest possible, for more complex operations, we would need to define them via `KernelAbstractions` kernels. We will not go into further details on that in this notebook.

However, now we have everything our model needs and we can finally use it! 

## Running our model 

First, we will define our initial conditions using `FieldInitializer`.

`FieldInitializer` can take in functions `(x,z)->val`, arrays or values. It uses `Oceananigans.set!(field, x)`, and allows all input arguments for `x` that `set!` allows: 
"""

# ╔═╡ fad517e2-9d3a-43cd-8e4d-0b72ca55b8c8
initializer = FieldInitializers(u = 1.0)

# ╔═╡ f2d02218-76f6-4b3a-84ca-38772f55d428
md"""
Then, we define our forcing. For that, our time-dependent forcing is loaded in from a `Oceananigans.FieldTimeSeries`. If you want to load the forcing from e.g. a netCDF file you can use the `RasterInputSource` that is based on `Rasters.jl`. In the concrete case, we'll just generate a random forcing: 
"""

# ╔═╡ 252af6a1-73c8-4abe-8100-690564641b0d
begin 
	t_F = 0:1:300
	F = FieldTimeSeries(grid, XY(), t_F)
	F.data .= randn(size(F));
	input = InputSource(; F)
end

# ╔═╡ 452f95e1-3c6b-4e49-935f-1a6f96c96bbb
md"""
Here we constructed a 2D (`XY()`) time series on our `grid` at times `t_F` with random normal distributed data and definted our `InputSource` for our model based on it. 
"""

# ╔═╡ 4b483c23-9e15-4d03-b275-8f530854669e
md"""
Then, we construct our model from the chosen `grid` and `initializer`
"""

# ╔═╡ 2a4234c5-f529-4166-94c3-0556565348ea
model = ExpModel(grid; initializer)

# ╔═╡ 4c36fdc0-5120-46b9-86ca-e875e23a6c1d
md"""
We now can initialize our model, i.e. we run all pre-computation, and initialize a numerical integrator for our model by passing it
to `initialize` along with a suitable timestepper and our input/forcing data, which we here choose to be the second-order Heun method with a timestep of 1 second.
"""

# ╔═╡ 7e38132b-d406-4863-b88f-90efe2a1bfa2
integrator = initialize(model, ForwardEuler(Δt = 1.0), input)

# ╔═╡ ab442662-9975-42e5-b5c7-48687f8cbe12
md"""
We can advance our model by one step via the `timestep!` method:
"""

# ╔═╡ 879d86d2-6828-4957-9aac-cd43508cbf1a
timestep!(integrator)

# ╔═╡ 4676ab3b-4f8f-4f47-9538-5f1e4ef257b1
integrator.state.u

# ╔═╡ 21e20c28-dfe1-4a0a-992f-c3499fbe4be8
md"""
or we can use `run!` for a fixed number of `steps` or over a desired `Dates.Period`:
"""

# ╔═╡ de3d4210-c39f-11f0-3d50-3f95a2361e2a
run!(integrator, period = Hour(1))

# ╔═╡ cce4d4d3-0fa4-4376-bcb6-c52603bc17d6
integrator.state.u

# ╔═╡ 7fa2dfbf-7077-4162-bcc1-ba2bd12b093c
md"""
As you can see, in just 1 hour of simulated time, our state variable already grew from `1` to `4e16`! If that's not exponential growth, we don't know what is ;) 
"""

# ╔═╡ 4c6d76e8-bc92-4abd-b2e8-15d26f5d4953
md"""
But wait there's more! What if we want to actually save the results?

The `integrator` data structure implements the Oceananigans model interface, so we can also use it to set up a `Simulation`:
"""

# ╔═╡ 95f479e2-2ffa-4e15-8952-421465eab2ee
sim = Simulation(integrator; stop_time = 300.0, Δt = 1.0)

# ╔═╡ 26000a4e-77cb-4c04-aeb2-ba5b0e14112a
begin
    # We need to import some types from Oceananigans here for output handling
    using Oceananigans: TimeInterval, JLD2Writer
    using Oceananigans.Units: seconds

    # Reset the integrator to its initial state
    Terrarium.initialize!(integrator)

    output_file = tempname()
    sim.output_writers[:snapshots] = JLD2Writer(
        integrator,
        (u = integrator.state.u,);
        filename = output_file,
        overwrite_existing = true,
        schedule = TimeInterval(10seconds)
    )
end

# ╔═╡ 081d0b29-927c-4a03-a3dd-4dcac043dcc1
md"""
We can then add an output writer to the simulation,
"""

# ╔═╡ 09118f2e-6c41-49e3-abf2-92b70976d755
md"""
and finally `run!` it!
"""

# ╔═╡ 4d416eb0-fbbc-4fec-939c-5c3909b2cef2
run!(sim)

# ╔═╡ 0f607788-53e7-4a55-95f0-3690e9867099
md"""
Then load the output data and plot the results:
"""

# ╔═╡ dbe8d0fa-893f-4c05-9e46-220ab41636f3
# Load output into field time series
fts = FieldTimeSeries(output_file, "u")

# ╔═╡ c06502ff-c021-488c-a333-36233091d046
plot(1:length(fts), [fts[i][1, 1, 1] for i in 1:length(fts)])

# ╔═╡ 25e22154-946f-4c32-a1fa-73d86e935ff3
md"""
Well that's it. We defined and ran a simple exponential model with external forcing following the Terrarium `AbstractModel` interface! Stay tuned for more!
"""

# ╔═╡ Cell order:
# ╟─5630efd5-2482-463d-913f-9addb120beec
# ╟─808d5d89-c1d2-4f6a-bd55-4b3a8444c90f
# ╠═94d82d31-42ec-41de-91e9-b5585b3a72d4
# ╟─07c8a3a4-21aa-4213-a876-eadc8754d4a0
# ╟─4922e264-c80d-4a5b-8891-a7c8a3fdbfe7
# ╠═78f268ef-5385-4c63-bc35-2c973de69da5
# ╟─054a8b11-250f-429f-966f-ca3c9a5dc2ef
# ╠═407786db-607f-4508-b697-fe75b3ce0b25
# ╟─575d920c-b12e-493f-95a7-5c962c3591fd
# ╠═82e45724-ba16-4806-9470-5cb4c43ea734
# ╟─d4d19de7-6f77-4873-9182-9832d1ca4381
# ╠═5ea313fc-3fbb-4092-a2cc-e0cd1f2fe641
# ╠═3815424f-6210-470d-aef1-99c60c71072f
# ╟─32373599-768f-4809-acdd-4704acc3f30b
# ╠═d55aaf4c-3033-45ba-9d64-8fa8ae4b671c
# ╠═5c8be7e4-f150-492b-a75d-96887a11f6da
# ╟─8d331856-6e9b-41d4-b1be-a84d5fedac8d
# ╠═fad517e2-9d3a-43cd-8e4d-0b72ca55b8c8
# ╟─f2d02218-76f6-4b3a-84ca-38772f55d428
# ╠═252af6a1-73c8-4abe-8100-690564641b0d
# ╟─452f95e1-3c6b-4e49-935f-1a6f96c96bbb
# ╟─4b483c23-9e15-4d03-b275-8f530854669e
# ╠═2a4234c5-f529-4166-94c3-0556565348ea
# ╟─4c36fdc0-5120-46b9-86ca-e875e23a6c1d
# ╠═7e38132b-d406-4863-b88f-90efe2a1bfa2
# ╟─ab442662-9975-42e5-b5c7-48687f8cbe12
# ╠═879d86d2-6828-4957-9aac-cd43508cbf1a
# ╠═4676ab3b-4f8f-4f47-9538-5f1e4ef257b1
# ╟─21e20c28-dfe1-4a0a-992f-c3499fbe4be8
# ╠═de3d4210-c39f-11f0-3d50-3f95a2361e2a
# ╠═cce4d4d3-0fa4-4376-bcb6-c52603bc17d6
# ╟─7fa2dfbf-7077-4162-bcc1-ba2bd12b093c
# ╟─4c6d76e8-bc92-4abd-b2e8-15d26f5d4953
# ╠═95f479e2-2ffa-4e15-8952-421465eab2ee
# ╟─081d0b29-927c-4a03-a3dd-4dcac043dcc1
# ╠═26000a4e-77cb-4c04-aeb2-ba5b0e14112a
# ╟─09118f2e-6c41-49e3-abf2-92b70976d755
# ╠═4d416eb0-fbbc-4fec-939c-5c3909b2cef2
# ╟─eeb283fa-5360-4bab-83cf-dcbc0bee7949
# ╟─0f607788-53e7-4a55-95f0-3690e9867099
# ╠═dbe8d0fa-893f-4c05-9e46-220ab41636f3
# ╠═c06502ff-c021-488c-a333-36233091d046
# ╟─25e22154-946f-4c32-a1fa-73d86e935ff3