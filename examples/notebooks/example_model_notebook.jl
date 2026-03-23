### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 808d5d89-c1d2-4f6a-bd55-4b3a8444c90f
begin
    import Pkg
    Pkg.activate(".")
    Pkg.develop(Pkg.PackageSpec(path = "../.."))
    Pkg.instantiate()
end

# ╔═╡ 94d82d31-42ec-41de-91e9-b5585b3a72d4
using Terrarium

# ╔═╡ 07c8a3a4-21aa-4213-a876-eadc8754d4a0
begin # for reproducability set a random seed
    using Random
    Random.seed!(1234)
end

# ╔═╡ eeb283fa-5360-4bab-83cf-dcbc0bee7949
using CairoMakie, GeoMakie

# ╔═╡ a55711ae-919c-46b0-a03e-ac4e105e0c4c
begin
    using RingGrids, NCDatasets
    global_grid = FullGaussianGrid(22) # we define the global grid we model on as a HealPIX grid
end

# ╔═╡ d77c8a4b-b53c-4906-a841-8ec37287ae9d
begin
    using PlutoUI
    PlutoUI.LocalResource("snow_depth.mp4")
end

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

The `grid` defines the spatial discretization. Our grids are based on those of [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) (and [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl)/[RingGrids.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl/tree/main/RingGrids)) in order to take advantage of their capabilities for device-agnostic parallelization.

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
# ╠═╡ disabled = true
#=╠═╡
grid = ColumnGrid(CPU(), Float64, UniformSpacing(N = 1))
  ╠═╡ =#

# ╔═╡ 054a8b11-250f-429f-966f-ca3c9a5dc2ef
md"""
## Defining the model

We start by defining a `struct` for our model that inherits from `AbstractModel` and consists of three properties: the spatial `grid`, an `initializer`, and a single `AbstractProcess` defining the dynamics, which we will also implement below.
"""

# ╔═╡ 407786db-607f-4508-b697-fe75b3ce0b25
@kwdef struct LinearDynamics{NF} <: Terrarium.AbstractProcess{NF}
    "Exponential growth rate"
    alpha::NF = 0.01

    "Constant offset for exponential growth"
    c::NF = 0.1
end

# ╔═╡ ecd92bff-a116-493d-9ce0-a2eb7d161dc6
@kwdef struct ExpModel{NF, Grid <: Terrarium.AbstractLandGrid{NF}, Dyn, Init} <: Terrarium.AbstractModel{NF, Grid}
    "Spatial grid on which state variables are discretized"
    grid::Grid

    "Linear dynamics resulting in exponential growth/decay"
    dynamics::Dyn = LinearDynamics()

    "Model initializer"
    initializer::Init = DefaultInitializer(eltype(grid))
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
    Terrarium.prognostic(:u, Terrarium.XY(), desc = "Exponential growth variable"),
    Terrarium.auxiliary(:c, Terrarium.XY(), desc = "Constant offset for exponential growth"),
    Terrarium.input(:F, Terrarium.XY(), desc = "External forcing"),
)

# ╔═╡ dfc52b4e-a015-4295-b47f-1dd2b10abeb2
begin
    using KernelAbstractions # we need this later for the kernels

    @kwdef struct SnowModel{NF, Grid <: Terrarium.AbstractLandGrid{NF}, Pro, Init} <: Terrarium.AbstractModel{NF, Grid}
        "Spatial grid on which state variables are discretized"
        grid::Grid

        "Snow melting process"
        snow_melt::Pro = DegreeDaySnow()

        "Model initializer"
        initializer::Init = DefaultInitializer(eltype(grid))
    end

    Terrarium.variables(model::SnowModel) = (
        Terrarium.variables(model.snow_melt)...,
    )

    @kwdef struct DegreeDaySnow{NF} <: Terrarium.AbstractProcess{NF}
        "Degree-day factor [m/(Ks)]"
        k::NF = 0.005f0 / (24 * 60 * 60)

        "Melting point of snow on the ground [°C]"
        T_melt::NF = 0.0f0
    end

    Terrarium.variables(model::DegreeDaySnow{NF}) where {NF} = (
        Terrarium.input(:air_temperature, XY(), default = NF(0), units = u"°C", desc = "Near-surface air temperature in °C"),
        Terrarium.input(:snow_fall, XY(), default = NF(0), units = u"m/s", desc = "snow fall rate in m/s"),
        Terrarium.prognostic(:snow_depth, XY(), units = u"m", desc = "Snow depth in m"),
    )

end

# ╔═╡ d4d19de7-6f77-4873-9182-9832d1ca4381
md"""
Here, we defined our three variables with their names as a `Symbol` and whether they are 2D variables (`XY`) on the spatial grid or 3D variables (`XYZ`) that also vary along the vertical z-axis. Here we are considering only a simple scalar model so we choose 2D (`XY`), bearing in mind that all points in the X and Y dimensions of `ColumnGrid` are independent of each other.

We also need to define `compute_auxiliary!` and `compute_tendencies!` as discussed above. We will use here a pattern which is commonly employed within Terrarium: we unpack the process from the model and forward the method calls to more specialzied ones defined for the `LinearDynamics` process. The `compute_auxiliary!` and `compute_tendencies!` of `AbstractProcess`es follow the signatures `(state, grid, processes...)`, as you see here: 
"""

# ╔═╡ 5ea313fc-3fbb-4092-a2cc-e0cd1f2fe641
function Terrarium.compute_auxiliary!(state, model::ExpModel)
    compute_auxiliary!(state, model.grid, model.dynamics)
    return nothing
end

# ╔═╡ 3815424f-6210-470d-aef1-99c60c71072f
function Terrarium.compute_tendencies!(state, model::ExpModel)
    compute_tendencies!(state, model.grid, model.dynamics)
    return nothing
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
        grid,
        dynamics::LinearDynamics
    )
    # set auxiliary variable for offset c
    return state.auxiliary.c .= dynamics.c
end

# ╔═╡ 5c8be7e4-f150-492b-a75d-96887a11f6da
# du/dt = u + c
function Terrarium.compute_tendencies!(
        state,
        grid,
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

First, we will define our initial conditions.

User-specified `Field` initializers passed to `initialize` can be provided in any form supported by `Oceananigans.set!`, including constants, arrays, and functions of the form `(x,z) -> val`: 
"""

# ╔═╡ fad517e2-9d3a-43cd-8e4d-0b72ca55b8c8
initializers = (u = 1.0,)

# ╔═╡ f2d02218-76f6-4b3a-84ca-38772f55d428
md"""
Then, we define our forcing. For that, our time-dependent forcing is loaded in from a `Oceananigans.FieldTimeSeries`. If you want to load the forcing from e.g. a netCDF file you can use the `RasterInputSource` that is based on [Rasters.jl](https://github.com/rafaqz/Rasters.jl). In the concrete case, we'll just generate a random forcing: 
"""

# ╔═╡ 252af6a1-73c8-4abe-8100-690564641b0d
#=╠═╡
begin
    t_F = 0:1:300
    F = FieldTimeSeries(grid, XY(), t_F)
    F.data .= randn(size(F))
    input = InputSource(grid, F, name=:F)
end
  ╠═╡ =#

# ╔═╡ 452f95e1-3c6b-4e49-935f-1a6f96c96bbb
md"""
Here we constructed a 2D (`XY()`) time series on our `grid` at times `t_F` with random normal distributed data and definted our `InputSource` for our model based on it. 
"""

# ╔═╡ 4b483c23-9e15-4d03-b275-8f530854669e
md"""
Then, we construct our model from the chosen `grid`
"""

# ╔═╡ 2a4234c5-f529-4166-94c3-0556565348ea
# ╠═╡ disabled = true
#=╠═╡
model = ExpModel(grid)
  ╠═╡ =#

# ╔═╡ 4c36fdc0-5120-46b9-86ca-e875e23a6c1d
md"""
We now can initialize our model, i.e. we run all pre-computation, and initialize a numerical integrator for our model by passing it
to `initialize` along with a suitable timestepper and our input/forcing data, which we here choose to be the second-order Heun method with a timestep of 1 second.
"""

# ╔═╡ 7e38132b-d406-4863-b88f-90efe2a1bfa2
# ╠═╡ disabled = true
#=╠═╡
integrator = initialize(model, Heun(Δt = 1.0), input; initializers)
  ╠═╡ =#

# ╔═╡ ab442662-9975-42e5-b5c7-48687f8cbe12
md"""
We can advance our model by one step via the `timestep!` method:
"""

# ╔═╡ 879d86d2-6828-4957-9aac-cd43508cbf1a
#=╠═╡
timestep!(integrator)
  ╠═╡ =#

# ╔═╡ 4676ab3b-4f8f-4f47-9538-5f1e4ef257b1
#=╠═╡
integrator.state.u
  ╠═╡ =#

# ╔═╡ 21e20c28-dfe1-4a0a-992f-c3499fbe4be8
md"""
or we can use `run!` for a fixed number of `steps` or over a desired `Dates.Period`:
"""

# ╔═╡ de3d4210-c39f-11f0-3d50-3f95a2361e2a
#=╠═╡
run!(integrator, period = Hour(1))
  ╠═╡ =#

# ╔═╡ cce4d4d3-0fa4-4376-bcb6-c52603bc17d6
#=╠═╡
integrator.state.u
  ╠═╡ =#

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
#=╠═╡
sim = Simulation(integrator; stop_time = 300.0, Δt = 1.0)
  ╠═╡ =#

# ╔═╡ 081d0b29-927c-4a03-a3dd-4dcac043dcc1
md"""
We can then add an output writer to the simulation and finally `run!` it!
"""

# ╔═╡ 26000a4e-77cb-4c04-aeb2-ba5b0e14112a
#=╠═╡
begin
    # We need to import some types from Oceananigans here for output handling
    using Oceananigans: TimeInterval, JLD2Writer
    using Oceananigans.Units: seconds

    # Reset the integrator to its initial state
    Terrarium.initialize!(integrator)

    output_dir = mkpath(tempname())
    output_file = joinpath(output_dir, "simulation.jld2")
    sim.output_writers[:snapshots] = JLD2Writer(
        integrator,
        (u = integrator.state.u,);
        filename = output_file,
        overwrite_existing = true,
        including = [:grid],
        schedule = TimeInterval(10seconds)
    )

    # Run the simulation
    run!(sim)
    @assert isfile(output_file) "Output file does not exist!"
    display("Simulaton data saved to $(output_file)")
end
  ╠═╡ =#

# ╔═╡ 0f607788-53e7-4a55-95f0-3690e9867099
md"""
Then load the output data and plot the results:
"""

# ╔═╡ 25e22154-946f-4c32-a1fa-73d86e935ff3
md"""
Well, this was a simple way to define and run a simple exponential model with external forcing following the Terrarium `AbstractModel` interface. 

But typically, our computations will be a bit more complicated than that and we can't easily assign them with a broadcastable operation like with done here in `compute_tendencies!`. So what do we have to do in these cases? 

## Writing kernelized-code for Terrarium 

Terrarium.jl is a device-agnostic modelling framework that runs across different architectures like x86 CPUs, ARM CPUs, but also GPUs. To achieve this wide compability, we rely on KernelAbstractions to write our compute code in kernels. But don't panic! We provide a lot of utilities and functions that help you with that and ensure that our models are fast and efficient as well. In simple terms, kernels can be thought of as the inner body of a for-loop. The kernel function implements one iteration of that loop, in Terrarium the kernel function implements the computation for a single column / grid point. To execute the computation the kernel is 'launched' on the device we set up when constructing the model (per default your CPU). In this launch we also set the range this kernel should iterate over (usually all columns of the grid), and hand over all arguments the kernel function needs. 

To demonstrate this process and all tools we have availble to help you with that, we will implement a degree day snow model. 

### Degree Day Model 

A degree day model models snow melt by assuming that snow depth decreases linearly  over time with the temperature when it is above its melting point. Following [Kavetski and Kuczera's formulation](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2006WR005195) we denote the snow mass balance as 

```math 
\frac{dS}{dt} = P - M
```

with snow storage $S \in \mathbb{R}^+$, the snow input $P$ and melt rate $M$ modelled as 

```math 
M = \begin{cases} 0 & \text{if } T \leq T_k, \\ 
	k(T - T_k) & \text{if } T > T_k \end{cases}
```

where $k$ is the degree-day factor/parameter and $T_k$ is a parameter for the melting point of snow on the ground. 

For our experiment we will use the degree day model with a prescribed surface temperature $T$, and initialize a global model with snow everywhere to watch the snow melt. 

### Your first kernel model 

First, we need to define our model that holds the snow melting process similar to before: 
"""

# ╔═╡ 52a2bf95-e258-41ab-922e-f0965d0d0ee2
md"""
Then, we need to define the dynamics again. Typically, we launch kernels on the level of `AbstractProcess`es in Terrarium. These processes then might have further parametrizations attached to them that need to be computed. We have a few utilities and conventions for this purpose:

* Kernels are typically named `compute_*_kernel!`, and follow a signature `compute_*_kernel!(output_fields, grid, fields, processes..., args...)`.
* In this function call we don't hand over the full `StateVariables` anymore as this would come with a large computational overhead on GPUs. Instead we only hand over those fields that the process actually needs. You can either assemble these fields manually, or use the convenience function `get_fields(state, processes...)` that returns all the fields that are defined in the `variables` of the respective processes or the model.
* Kernels are launched with the `launch!` function that defines whether the kernel is launched over all columns / gridpoints `XY` or over all columns and horizontal layers `XYZ` of a `grid`
* Kernel functions themselves are supposed to be very minimal functions that forward the actual computation to pointwise compute functions that compute the needed quantity for the grid points `i,j`. This increases the reusability and composability of our model code. These compute functions are really the core of the implementation of our model and similar to those you typically find in more traditional land models: They compute the action of a process for a single grid point given inputs and parameters. 
* The compute functions follow either a mutating pattern `compute_*!(out, i, j, grid, fields, processes..., args...)` or a non-mutating pattern `compute_*(i, j, grid, fields, processes..., args...)`

Let's put all of this in practice now for our example. 

First, we need to define the `compute_tendencies!` for our `Model` as before, then we launch the kernel in the `compute_tendencies!` of our `Process`. For this model we don't actually have auxiliary variables, so we don't have to actually define a `compute_auxiliary!` in this case, as we inherit a default `compute_auxiliary!(args...) = nothing`. 
"""

# ╔═╡ ab4a216f-3962-4a6f-8f92-2e9a08798e7c
md"""
The `@index(Global, NTuple)` is a function from KernelAbstractions that gets us the current index of our iteration, so the column we compute. 

With that, we can also define our `compute_foo_tendency` function that does the actual computation: 
"""

# ╔═╡ d8b05ae3-ecba-41de-84c7-45cbf31b735d
function compute_snow_melt_tendency(i, j, grid, fields, snow_melt)
    # get the variables we need
    P = fields.snow_fall[i, j]
    T = fields.air_temperature[i, j]

    # get the parameters
    T_melt = snow_melt.T_melt
    k = snow_melt.k

    if T > T_melt
        return P - k * (T - T_melt)
    else
        return P
    end
end

# ╔═╡ 72ec62c7-5066-481d-b9c9-84a4851a1e0c
begin
    function Terrarium.compute_tendencies!(state, model::SnowModel)
        compute_tendencies!(state, model.grid, model.snow_melt)
        return nothing
    end

    function Terrarium.compute_tendencies!(
            state,
            grid,
            snow_melt::DegreeDaySnow
        )
        fields = get_fields(state, snow_melt)
        return Terrarium.launch!(grid, XY, compute_snow_melt!, state.tendencies, fields, snow_melt)
    end

    @kernel function compute_snow_melt!(tend, grid, fields, snow_melt)
        i, j = @index(Global, NTuple)
        tend.snow_depth[i, j] = compute_snow_melt_tendency(i, j, grid, fields, snow_melt)
    end

    # no auxiliary variables
    Terrarium.compute_auxiliary!(state, model::SnowModel) = nothing
    Terrarium.compute_auxiliary!(state, grid, model::DegreeDaySnow) = nothing
end

# ╔═╡ e81f4b38-5789-416e-acee-e02b052cb8f4
md"""
For a simple process like this, this is of course quite a lot of overhead, but this structure allows us to also compose very complex land models with ease and efficiency. 

There's one additional thing we need to take of: the snow depth is strictly non-negative, but a simple model like ours implemented would quickly reach negative values for $S$. While the aforementioned [Kavetski and Kuczera paper](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2006WR005195) mitigates this issue by smoothing the model equation, we will in this example do the simplest possible strategy: we simply clip non-negative values during the time stepping of our model. To implement such a clipping we have to extend `timestep!(state::StateVaribales, model::ExpModel, timestepper, Δt)`. This method is applied after each explicit step the time stepper takes (but before any closure relations are applied if they exist). Let's implement the clipping: 
"""

# ╔═╡ b723c568-c0e1-4d9a-9a74-237d7cfd1ea9
function Terrarium.timestep!(state::StateVariables, model::ExpModel, timestepper, Δt)
    return state.snow_depth .= max.(state.snow_depth, 0)
end

# ╔═╡ 841c540f-ed63-4d89-9baf-836ccb3aed0d
md"""
But, now we really have implemented everything we need for our dynamics and we can finally run the model again. For this we set up our initializers again in a very similar manner as above for the previous example, just this time with inputs from netCDF files and a global grid: 
"""

# ╔═╡ 3ef9f3c7-16a2-416f-981a-64427d89b033
md"""
Now, we get all the input data. We can use the input data provided by `RingGrids` / `SpeedyWeather` in this case: 
"""

# ╔═╡ 86f4103b-ddaf-4f89-93cf-1290c623274e
begin
    snow_climatology = RingGrids.get_asset("data/boundary_conditions/snow.nc", from_assets = true, name = "snow", ArrayType = FullGaussianField, FileFormat = NCDataset, output_grid = global_grid) ./ 3.8e10 # ~ conversion from kg/(month * m^2) to m/(s * m^2)

    lst_climatology = RingGrids.get_asset("data/boundary_conditions/land_surface_temperature.nc", from_assets = true, name = "lst", ArrayType = FullGaussianField, FileFormat = NCDataset, output_grid = global_grid) .- 273.15 # data is in K, we want C

    # the land sea mask we infer from the non-NaN points in the climatology files
    land_sea_mask = isfinite.(snow_climatology[:, 1])
    @assert all(land_sea_mask .== isfinite.(lst_climatology[:, 1])) # make sure it's the same for both

end

# ╔═╡ 6ba13cea-da1c-457e-ada0-8987a0667b24
md"""
The snow and land surface temperatures are monthly climatologies. For this simple example, we'll just pick the January value. Let's quickly look at our input data. 
"""

# ╔═╡ fc8562fd-0213-48be-870f-ab9b06c54543
heatmap(land_sea_mask)

# ╔═╡ 051d99da-a7b8-4cd8-be7c-3f7f615345d3
heatmap(lst_climatology[:, 1], title = "Land Surface Temperature")

# ╔═╡ 9982a3fd-c2ef-4bba-917e-211912fdce84
heatmap(snow_climatology[:, 1], title = "Snow")

# ╔═╡ 49364e74-1272-4d65-892c-5b08d0e49a54
md"""
Ok, so now let's put everything together! 

* We defined our model `SnowModel` and dynamics `DegreeDaySnow`
* We loaded climatological input data and a land sea mask for our grid 

Now, we just need to define initialize everything correctly. As we are working with globally gridded data, we will define `ColumnRingGrid` based on the `global_grid` we already initialized. Then, we will load our inputs. For this, we will choose the January (so the first element) of our climatology files. When using them in `InputSource` be sure to choose the same name and units as used in the definitions of the dynamics before.
"""

# ╔═╡ e80009fb-cf05-4360-9f7e-c355d059ff5c
begin
    snow_grid = ColumnRingGrid(UniformSpacing(N = 1), global_grid, land_sea_mask)
    snow_input = InputSource(snow_grid, snow_climatology[:, 1], name = :snow_fall, units = u"m/s")
    lst_input = InputSource(snow_grid, lst_climatology[:, 1], name = :air_temperature, units = u"°C")
end

# ╔═╡ 7f310793-f4af-4013-b02f-8273107246e7
md"""
As an initial condition, we just cover the whole Earth in deep snow (everywhere the same)!
"""

# ╔═╡ b12f3815-89d3-4c1d-9224-6bde2d1f939e
snow_initializers = (snow_depth = 0.5,)

# ╔═╡ 01e757f2-55ec-494f-8286-29e5e5fe0a2e
md"""
Now, we initialize our model and the integrator. As in the first example, we use a `Heun` time stepper
"""

# ╔═╡ a29583b9-62be-42ef-adf0-867a734a03d7
snow_model = SnowModel(snow_grid)

# ╔═╡ 018375eb-fd00-48b6-9b9d-dc42ba2d5c2d
snow_integrator = initialize(snow_model, Heun(Δt = Float32(1.0)), snow_input, lst_input; initializers = snow_initializers)

# ╔═╡ 5f9a57c7-9094-418f-ae62-11cce5cad690
md"""
... and we can finally run the model. As before, by wrapping it in an `Oceananigans.Simulation` to output our results 
"""

# ╔═╡ 4288f181-fb84-4ce1-b13b-674a5a8de132
snow_sim = Simulation(snow_integrator; stop_time = 7.0e6, Δt = 3600.0)

# ╔═╡ 7094e65a-fbd5-4bc1-9734-54c77fbeb757
begin
    # We need to import some types from Oceananigans here for output handling
    using Oceananigans: TimeInterval, JLD2Writer
    using Oceananigans.Units: seconds

    # Reset the integrator to its initial state
    Terrarium.initialize!(snow_integrator)

    output_dir = mkpath(tempname())
    output_file = joinpath(output_dir, "ddsnow-simulation.jld2")
    snow_sim.output_writers[:snapshots] = JLD2Writer(
        snow_integrator,
        (snow_depth = snow_integrator.state.snow_depth,);
        filename = output_file,
        overwrite_existing = true,
        including = [:grid],
        schedule = TimeInterval(3600) # output every hour
    )
end

# ╔═╡ dbe8d0fa-893f-4c05-9e46-220ab41636f3
# Load output into field time series
fts = FieldTimeSeries(output_file, "u")

# ╔═╡ c06502ff-c021-488c-a333-36233091d046
plot(1:length(fts), [fts[i][1, 1, 1] for i in 1:length(fts)])

# ╔═╡ bc5c7603-314a-411e-8097-a6344f7bf52a
begin
    using JLD2

    output_file_name = joinpath(output_dir, "ddsnow-simulation.jld2")
    fts_result = FieldTimeSeries(output_file, "snow_depth")
end

# ╔═╡ c8a89e9e-d24a-415c-9a5b-5089053f6384
begin
    # Run the simulation
    run!(snow_sim)
    @assert isfile(output_file) "Output file does not exist!"
    display("Simulaton data saved to $(output_file)")
end

# ╔═╡ b7c37a45-b00f-4d27-bcf5-f42ac610566e
md"""
And now, we plot that data. First we load the JLD2 file. 
"""

# ╔═╡ 78b731fd-6106-4a43-a6f5-264f5fbc271d
md""" 
Then, we plot it using `CairoMakie`. For this purpose we first convert to a `RingGrids.Field` and then plot it via `heatmap`
"""

# ╔═╡ 0872d1a8-f170-46ee-a99e-fe1c6f533ace
begin
    tsteps = 1

    ring_field = RingGrids.Field(fts_result[tsteps], snow_grid)[:, 1]
    heatmap(ring_field)

    fig = Figure(size = (1200, 660))

    ax = Axis(
        fig[1, 1],
        aspect = 2,             # 0-360˚E -90-90˚N maps have an aspect of 2:1
        title = "Snow Depth",
        titlesize = 10,
        xticks = 0:60:360,      # label 0˚E, 60˚E, 120˚E, ...
        yticks = -60:30:60,     # label -60˚N, -30˚N, 0˚N, ...
        xticklabelsize = 10,
        yticklabelsize = 10,
        xtickformat = values -> ["$(round(Int, value))˚E" for value in values],
        ytickformat = values -> ["$(round(Int, value))˚N" for value in values],
    )

    lond = RingGrids.get_lond(ring_field)    # get lon, lat axes in degrees
    latd = RingGrids.get_latd(ring_field)

    n_t = Observable(1)

    data = @lift Matrix(RingGrids.Field(fts_result[$n_t], snow_grid)[:, 1])
    hm = heatmap!(ax, lond, latd, data, colorrange = (0, 1))
    Colorbar(fig[:, end + 1], hm)

    frames = 1:length(fts_result)

    record(fig, "snow_depth.mp4", frames, framerate = 12) do i # core animation loop
        n_t[] = i
    end
end

# ╔═╡ bfb60a6c-df9c-4d45-88f0-a01f572fe8b2
md"""
And just like that we have implemented our snow simulation. In this version the forcing / input is completly static, so we converge to a static point that corresponds to the January climatology. That's why still see a faily big snow cover in the northern hemisphere. An obvious next step for this model would be now to actually use the full seasonal climatology.
"""

# ╔═╡ Cell order:
# ╟─5630efd5-2482-463d-913f-9addb120beec
# ╠═808d5d89-c1d2-4f6a-bd55-4b3a8444c90f
# ╠═94d82d31-42ec-41de-91e9-b5585b3a72d4
# ╠═07c8a3a4-21aa-4213-a876-eadc8754d4a0
# ╟─4922e264-c80d-4a5b-8891-a7c8a3fdbfe7
# ╠═78f268ef-5385-4c63-bc35-2c973de69da5
# ╟─054a8b11-250f-429f-966f-ca3c9a5dc2ef
# ╠═407786db-607f-4508-b697-fe75b3ce0b25
# ╠═ecd92bff-a116-493d-9ce0-a2eb7d161dc6
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
# ╟─eeb283fa-5360-4bab-83cf-dcbc0bee7949
# ╟─0f607788-53e7-4a55-95f0-3690e9867099
# ╠═dbe8d0fa-893f-4c05-9e46-220ab41636f3
# ╠═c06502ff-c021-488c-a333-36233091d046
# ╟─25e22154-946f-4c32-a1fa-73d86e935ff3
# ╠═dfc52b4e-a015-4295-b47f-1dd2b10abeb2
# ╟─52a2bf95-e258-41ab-922e-f0965d0d0ee2
# ╠═72ec62c7-5066-481d-b9c9-84a4851a1e0c
# ╟─ab4a216f-3962-4a6f-8f92-2e9a08798e7c
# ╠═d8b05ae3-ecba-41de-84c7-45cbf31b735d
# ╟─e81f4b38-5789-416e-acee-e02b052cb8f4
# ╠═b723c568-c0e1-4d9a-9a74-237d7cfd1ea9
# ╟─841c540f-ed63-4d89-9baf-836ccb3aed0d
# ╠═a55711ae-919c-46b0-a03e-ac4e105e0c4c
# ╟─3ef9f3c7-16a2-416f-981a-64427d89b033
# ╠═86f4103b-ddaf-4f89-93cf-1290c623274e
# ╟─6ba13cea-da1c-457e-ada0-8987a0667b24
# ╠═fc8562fd-0213-48be-870f-ab9b06c54543
# ╠═051d99da-a7b8-4cd8-be7c-3f7f615345d3
# ╠═9982a3fd-c2ef-4bba-917e-211912fdce84
# ╟─49364e74-1272-4d65-892c-5b08d0e49a54
# ╠═e80009fb-cf05-4360-9f7e-c355d059ff5c
# ╟─7f310793-f4af-4013-b02f-8273107246e7
# ╠═b12f3815-89d3-4c1d-9224-6bde2d1f939e
# ╠═01e757f2-55ec-494f-8286-29e5e5fe0a2e
# ╠═a29583b9-62be-42ef-adf0-867a734a03d7
# ╠═018375eb-fd00-48b6-9b9d-dc42ba2d5c2d
# ╟─5f9a57c7-9094-418f-ae62-11cce5cad690
# ╠═4288f181-fb84-4ce1-b13b-674a5a8de132
# ╠═7094e65a-fbd5-4bc1-9734-54c77fbeb757
# ╠═c8a89e9e-d24a-415c-9a5b-5089053f6384
# ╟─b7c37a45-b00f-4d27-bcf5-f42ac610566e
# ╠═bc5c7603-314a-411e-8097-a6344f7bf52a
# ╟─78b731fd-6106-4a43-a6f5-264f5fbc271d
# ╠═0872d1a8-f170-46ee-a99e-fe1c6f533ace
# ╠═d77c8a4b-b53c-4906-a841-8ec37287ae9d
# ╟─bfb60a6c-df9c-4d45-88f0-a01f572fe8b2
