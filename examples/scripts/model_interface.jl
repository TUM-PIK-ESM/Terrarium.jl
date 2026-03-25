# # Building a Land Model with Terrarium
#
# In this example we demonstrate Terrarium's model interface by implementing two models:
# a simple exponential growth model (introducing the interface) and a global degree day
# snow model that uses GPU-compatible kernels.

using Terrarium
using CairoMakie, GeoMakie
using KernelAbstractions
using RingGrids, NCDatasets
using Oceananigans: TimeInterval, JLD2Writer
using Oceananigans.Units: seconds
using JLD2
using Random
Random.seed!(1234)

# ## A Simple Exponential Growth Model
#
# In this example we will set up an embarrassingly simple example to demonstrate Terrarium's model interface. Our model will have 1-dimensional exponential dynamics with a constant offset
#
#```math
#\frac{du}{dt} = \alpha u + c + F(t)
#```
#
# for an arbitrary prognostic variable ``u``. For the sake of this demonstration we will treat the offset ``c`` as an auxiliary/diagnostic variable even though it is constant in time. ``F(t)`` is an external forcing that we apply.
#
# We begin by defining our model `struct` that subtypes `Terrarium.AbstractModel`:
#
# A "model" in Terrarium is a subtype of `Terrarium.AbstractModel` and is a `struct` type constisting of
#  * `grid` which defines the discretization of the spatial domain
#  * `initializer` which is responsible for initializing state variables
#  * further fields that define processes, dynamics and submodels
#
# When we follow the advised naming notations of `grid` and `initializer` we inherit default methods from `Terrarium.AbstractModel` such as `get_grid` and `get_initializer`. For more complex models we might need to implement custom overrides of `initialize!(state, ::Model, ::Initializer)` to initialize model states.
#
# ## What is a "grid"?
#
# The `grid` defines the spatial discretization. Our grids are based on those of [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) (and [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl)/[RingGrids.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl/tree/main/RingGrids)) in order to take advantage of their capabilities for device-agnostic parallelization.
#
# Terrarium currently provides two grid types:
#
# * `ColumnGrid` is a set of laterally independent vertical columns with dimensions ``(x, y, z)`` where ``x`` is the column dimension, ``y=1`` is constant, and ``z`` is the vertical axis,
# * `ColumnRingGrid` represents a global (spherical) grid of independent, vertical columns where the spatial discretization in the horizontal direction is defined by a `RingGrids.AbstractGrid`.
#
# In both cases we need to specificy the vertical discretizataion via an `UniformSpacing`, `ExponentialSpacing` or `PrescribedSpacing`.
#
# ## Initializer and Boundary Conditions
#
# For our basic example here the default initializer (which does nothing) will suffice, and we won't have to define a custom one.
#
# Boundary conditions are specified by passing Oceananigans `BoundaryCondition` types to `initialize`. In the case of a linear ODE, however, no boundary conditions are required.
#
# ## What's our `grid`?
#
# For our current example, we are defining a simple linear ODE without any spatial dynamics, so we can get away with just a single column with one vertical layer. We can define it like so:

grid = ColumnGrid(Terrarium.CPU(), Float64, UniformSpacing(N = 1))

# ## Defining the model
#
# We start by defining a `struct` for our model that inherits from `AbstractModel` and consists of three properties: the spatial `grid`, an `initializer`, and a single `AbstractProcess` defining the dynamics, which we will also implement below.

@kwdef struct LinearDynamics{NF} <: Terrarium.AbstractProcess{NF}
    "Exponential growth rate"
    alpha::NF = 0.01
    "Constant offset"
    c::NF = 0.1
end
#
@kwdef struct ExpModel{NF, Grid <: Terrarium.AbstractLandGrid{NF}, Dyn, Init} <: Terrarium.AbstractModel{NF, Grid}
    "Spatial grid on which state variables are discretized"
    grid::Grid
    "Linear dynamics process"
    dynamics::Dyn = LinearDynamics()
    "Model initializer"
    initializer::Init = DefaultInitializer(eltype(grid))
end

# ## Defining the model behaviour
#
# Now, we want to define our intended model behaviour. For this, we need to define the following methods:
#
# * `variables(::Model)` returns a tuple of variable metadata declaring the state variables. Variables must be one of three types: `prognostic`, `auxiliary` (sometimes referred to as "diagnostic"), or `input`. Prognostic variables fully characterize the state of the system at any given timestep and are updated according to their tendencies (i.e. ``u`` in our example). Tendencies are automatically allocated for each prognostic variable declared by the model. In this example we will treat the offset ``c`` as an auxiliary variable, though we could also just include it as a constant in the tendency computations.
# * `compute_auxiliary!(state, ::Model)` computes the values of all auxiliary variables (if necessary) assuming that the prognostic variables of the system in state are available for the current timestep.
# * `compute_tendencies!(state, ::Model)` computes the tendencies based on the current values of the prognostic and auxiliary variables stored in state.
#
# So, let's define those:

Terrarium.variables(::ExpModel) = (
    Terrarium.prognostic(:u, Terrarium.XY(), desc = "Exponential growth variable"),
    Terrarium.auxiliary(:c, Terrarium.XY(), desc = "Constant offset for growth"),
    Terrarium.input(:F, Terrarium.XY(), default = 0.0, desc = "External forcing"),
)

# Here, we defined our three variables with their names as a `Symbol` and whether they are 2D variables (`XY`) on the spatial grid or 3D variables (`XYZ`) that also vary along the vertical z-axis. Here we are considering only a simple scalar model so we choose 2D (`XY`), bearing in mind that all points in the X and Y dimensions of `ColumnGrid` are independent of each other.
#
# We also need to define `compute_auxiliary!` and `compute_tendencies!` as discussed above. We will use here a pattern which is commonly employed within Terrarium: we unpack the grid and process from the model and forward the method calls to more specialized ones defined for the `LinearDynamics` process. The `compute_auxiliary!` and `compute_tendencies!` of `AbstractProcess`es follow the signatures `(state, grid, processes...)`, as you see here:

function Terrarium.compute_auxiliary!(state, model::ExpModel)
    compute_auxiliary!(state, model.grid, model.dynamics)
    return nothing
end
#
function Terrarium.compute_tendencies!(state, model::ExpModel)
    compute_tendencies!(state, model.grid, model.dynamics)
    return nothing
end

# Note that, when implementing models within the Terrarium module itself, the `Terrarium.` qualifier in the definition is not needed.
#
# ## Implementing the dynamics
#
# Next, we define the functions that compute the actual dynamics. In order to do this, we need to know a little about how the variables we just defined are handled in our `StateVariables`. The `StateVariables` hold all prognostic and auxiliary variables, their tendencies and closures and additional inputs and forcings in seperate `NamedTuples`. Note that Terrarium also defines shortcuts such that, e.g. in our example, both `state.prognostic.u` and `state.u` would work.
#
# With that in mind, let's define the methods:

function Terrarium.compute_auxiliary!(state, grid, dynamics::LinearDynamics)
    return state.auxiliary.c .= dynamics.c
end
#
function Terrarium.compute_tendencies!(state, grid, dynamics::LinearDynamics)
    return let u = state.prognostic.u,
            ∂u∂t = state.tendencies.u,
            α = dynamics.alpha,
            c = state.auxiliary.c,
            F = state.inputs.F
        ∂u∂t .= α .* u .+ c .+ F
    end
end

# These example compute functions are really the simplest possible, for more complex operations, we would need to define them via `KernelAbstractions` kernels. We will not go into further details on that in this notebook.
#
# However, now we have everything our model needs and we can finally use it!
#
# ## Running our model
#
# First, we will define our initial conditions.
#
# User-specified `Field` initializers passed to `initialize` can be provided in any form supported by `Oceananigans.set!`, including constants, arrays, and functions of the form `(x,z) -> val`:

initializers = (u = 1.0,)

# Then, we define our forcing. For that, our time-dependent forcing is loaded in from a `Oceananigans.FieldTimeSeries`. If you want to load the forcing from e.g. a netCDF file you can use the `RasterInputSource` that is based on [Rasters.jl](https://github.com/rafaqz/Rasters.jl). In the concrete case, we'll just generate a random forcing:

t_F = 0:1:300;
F = FieldTimeSeries(grid, XY(), t_F);
F.data .= randn(size(F));
input = InputSource(grid, F, name = :F)

# Here we constructed a 2D (`XY()`) time series on our `grid` at times `t_F` with random normal distributed data and definted our `InputSource` for our model based on it.

# Then, we construct our model from the chosen `grid`

model = ExpModel(grid)

# We now can initialize our model, i.e. we run all pre-computation, and initialize a numerical integrator for our model by passing it
# to `initialize` along with a suitable timestepper and our input/forcing data, which we here choose to be the second-order Heun method with a timestep of 1 second.

integrator = initialize(model, Heun(Δt = 1.0), input; initializers)

# We can advance our model by one step via the `timestep!` method:

timestep!(integrator)

integrator.state.u

# As you can see, in just 1 hour of simulated time, our state variable already grew from `1` to `4e16`! If that's not exponential growth, we don't know what is ;)

# But wait there's more! What if we want to actually save the results?
#
# The `integrator` data structure implements the Oceananigans model interface, so we can also use it to set up a `Simulation`:

sim = Simulation(integrator; stop_time = 300.0, Δt = 1.0)

# We can then add an output writer to the simulation and finally `run!` it!

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
run!(sim)
@assert isfile(output_file) "Output file does not exist!"
println("Simulation data saved to $(output_file)")

# Then load the output data and plot the results:

fts = FieldTimeSeries(output_file, "u")

plot(1:length(fts), [fts[i][1, 1, 1] for i in 1:length(fts)])

# We have seen a simple example of how to define and run an exponential growth model with external forcing following the Terrarium `AbstractModel` interface.
#
# But typically, our computations will be a bit more complicated than that, and we can't just rely on simple broadcasting operations like what we did in `compute_tendencies!` above. The reason for this is simply efficiency: it is (usually) more efficient to bundle together many scalar operations into operations that can be massively parallelized. But how do we actually achieve this?
#
# ## Writing kernelized-code for Terrarium
#
# Terrarium.jl is a device-agnostic modelling framework that runs across different architectures like x86 CPUs, ARM CPUs, but also most importantly, GPUs. To achieve this wide compatibility, we rely on KernelAbstractions to turn our heavy computations into parallelizable kernels. But don't panic! We provide a lot of utilities and functions that help make this easy and ensure that the resulting models are fast and efficient as well.
#
# In simple terms, kernels can be thought of as the inner body of a for-loop. The kernel function implements one iteration of that loop; in Terrarium, the kernel function implements the computation for a single column / grid point. To execute the computation the kernel is "launched" on the device we set up when constructing the model (per default your CPU). When configuring the kernel launch, we also set the range this kernel should iterate over (usually all points on the grid), and hand over all arguments required by the kernel function.
#
# To demonstrate this process, and all tools we have available to facilitate it, we will implement a simplified degree day snow model.
#
# ### Degree Day Model
#
# A degree day model models snow melt by assuming that snow storage or snow water equivalent decreases linearly  over time with the temperature when it is above its melting point. Following [Kavetski and Kuczera's formulation](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2006WR005195) we denote the snow mass balance as
#
# ```math
# \frac{dS}{dt} = P - M
# ```
#
# with snow storage ``S \in \mathbb{R}^+``, the snow input ``P`` and melt rate ``M`` modelled as
#
# ```math
# M = \begin{cases} 0 & \text{if } T \leq T_k, \\
# 	k(T - T_k) & \text{if } T > T_k \end{cases}
# ```
#
# where ``k`` is the degree-day factor/parameter and ``T_k`` is a parameter for the melting point of snow on the ground.
#
# For our experiment we will use the degree day model with a prescribed surface temperature ``T``, and initialize a global model with snow everywhere to watch the snow melt.
#
# ### Your first kernel model
#
# First, we need to define our model that holds the snow melting process similar to how we constructed the `ExpModel`:

@kwdef struct DegreeDaySnow{NF} <: Terrarium.AbstractProcess{NF}
    "Degree-day factor [m/(Ks)]"
    k::NF = 0.005f0 / (24 * 60 * 60)
    "Melting point of snow on the ground [°C]"
    T_melt::NF = 0.0f0
end
#
Terrarium.variables(model::DegreeDaySnow{NF}) where {NF} = (
    Terrarium.input(:air_temperature, XY(), default = NF(0), units = u"°C", desc = "Near-surface air temperature in °C"),
    Terrarium.input(:snow_fall, XY(), default = NF(0), units = u"m/s", desc = "snow fall rate in m/s"),
    Terrarium.prognostic(:snow_storage, XY(), units = u"m", desc = "Snow water equivalent in m"),
)

@kwdef struct SnowModel{NF, Grid <: Terrarium.AbstractLandGrid{NF}, Pro, Init} <: Terrarium.AbstractModel{NF, Grid}
    "Spatial grid on which state variables are discretized"
    grid::Grid
    "Snow melting process"
    snow_melt::Pro = DegreeDaySnow()
    "Model initializer"
    initializer::Init = DefaultInitializer(eltype(grid))
end
#
Terrarium.variables(model::SnowModel) = (
    Terrarium.variables(model.snow_melt)...,
)

# Then, we need to define the dynamics. Typically, we launch kernels on the level of `AbstractProcess`es in Terrarium. These processes then might have further parameterizations attached to them that need to be computed. We have a few utilities and conventions for this purpose:
#
# * Kernels are typically named `compute_*_kernel!`, and follow a signature `compute_*_kernel!(output_fields, grid, fields, processes..., args...)`.
# We always hand over only the minimum set of `fields` that the process (and its dependencies) actually need. You can either assemble these fields manually, or use the convenience function `get_fields(state, processes...)` that returns all the fields that are defined in the `variables` of the respective processes or the model. Only the minimum set of `fields` is used, as handing over the full `state` to the kernels would come with high computational overhead due to the need to copy type information from the host device to the GPU.
# * Kernels are launched with the `launch!` function that defines whether the kernel is launched over all columns/grid points `XY` (2D) or also over all vertical layers `XYZ` (3D) of the model `grid`.
# * The functions that define the kernels with `@kernel` are supposed to be very minimal functions that forward the actual computation to pointwise compute functions that compute the needed quantity for the grid points `i,j`. This increases the reusability and composability of our model code. These compute functions are really the core of the implementation of our model and similar to functions you might find in more traditional land model, called from within a loop.They compute the action of a process for a single grid point given inputs and parameters.
# * The compute functions follow either a mutating pattern `compute_*!(out, i, j, grid, fields, processes..., args...)` or a non-mutating pattern `compute_*(i, j, grid, fields, processes..., args...)`. As a general rule, the mutating functions should defer all actual computations to non-mutating functions and then store the results in the `Field`s given in the first argument.
#
# Let's put all of this in practice now for our example.
#
# First, we need to define the `compute_tendencies!` for our `Model` as before, then we launch the kernel in the `compute_tendencies!` of our `Process`. For this model we don't actually have auxiliary variables, so we don't have to actually define a `compute_auxiliary!` in this case, as we inherit a default `compute_auxiliary!(args...) = nothing`.

function Terrarium.compute_tendencies!(state, model::SnowModel)
    compute_tendencies!(state, model.grid, model.snow_melt)
    return nothing
end
#
function Terrarium.compute_tendencies!(
        state,
        grid,
        snow_melt::DegreeDaySnow
    )
    fields = get_fields(state, snow_melt)
    return Terrarium.launch!(grid, XY, compute_snow_flux!, state.tendencies, fields, snow_melt)
end
#
@kernel function compute_snow_flux!(tend, grid, fields, snow_melt)
    i, j = @index(Global, NTuple)
    tend.snow_storage[i, j] = compute_snow_flux_tendency(i, j, grid, fields, snow_melt)
end
#
Terrarium.compute_auxiliary!(state, model::SnowModel) = nothing
Terrarium.compute_auxiliary!(state, grid, model::DegreeDaySnow) = nothing

# The `@index(Global, NTuple)` is a function from KernelAbstractions that gets us the current index of our iteration, so the column we compute.
#
# With that, we can also define our `compute_foo_tendency` function that does the actual computation:

function compute_snow_flux_tendency(i, j, grid, fields, snow_melt)
    ## get the variables we need
    P = fields.snow_fall[i, j]
    T = fields.air_temperature[i, j]
    ## get the parameters
    T_melt = snow_melt.T_melt
    k = snow_melt.k
    return ifelse(T > T_melt, P - k * (T - T_melt), P)
end

# For a simple process like this, this is of course quite a lot of overhead, but this structure allows us to very efficiently build up complex land models from relatively simple components.
#
# There's one additional thing we need to take care of: the snow storage is strictly non-negative, but our model as currently implemented would quickly reach negative values for ``S``. While the aforementioned [Kavetski and Kuczera paper](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2006WR005195) mitigates this issue by smoothing the model equation, we will in this example do the simplest possible strategy: we simply clip non-negative values during the time stepping of our model. To implement such a clipping we have to extend `timestep!(state::StateVaribales, model::ExpModel, timestepper, Δt)`. This method is applied after each explicit step the time stepper takes (but before any closure relations are applied if they exist). Let's implement the clipping:

function Terrarium.timestep!(state::StateVariables, model::SnowModel, timestepper::Terrarium.AbstractTimeStepper, Δt)
    interior(state.snow_storage) .= max.(interior(state.snow_storage), 0)
    return nothing
end

# But, now we really have implemented everything we need for our dynamics and we can finally run the model again. For this we set up our initializers again in a very similar manner as above for the previous example, just this time with inputs from netCDF files and a global grid:

global_grid = FullGaussianGrid(22) # we define the global grid we model on as a HealPIX grid

# Now, we get all the input data. We can use the input data provided by `RingGrids` / `SpeedyWeather` in this case:

snow_climatology = RingGrids.get_asset("data/boundary_conditions/snow.nc", from_assets = true, name = "snow", ArrayType = FullGaussianField, FileFormat = NCDataset, output_grid = global_grid) ./ 3.8e10; # ~ conversion from kg/(month * m^2) to m/(s * m^2)
lst_climatology = RingGrids.get_asset("data/boundary_conditions/land_surface_temperature.nc", from_assets = true, name = "lst", ArrayType = FullGaussianField, FileFormat = NCDataset, output_grid = global_grid) .- 273.15; # data is in K, we want C
land_sea_mask = isfinite.(snow_climatology[:, 1])
@assert all(land_sea_mask .== isfinite.(lst_climatology[:, 1]))

# The snow and land surface temperatures are monthly climatologies. For this simple example, we'll just pick the January value. Let's quickly look at our input data. First the land sea mask:

heatmap(land_sea_mask, title = "Land Sea Mask")

# Then, the land surface temperature:

heatmap(lst_climatology[:, 1], title = "Land Surface Temperature")

# And finally, the snowfall:

heatmap(snow_climatology[:, 1], title = "Snowfall")

# Ok, so now let's put everything together!
#
# * We defined our model `SnowModel` and dynamics `DegreeDaySnow`
# * We loaded climatological input data and a land sea mask for our grid
#
# Now, we just need to define initialize everything correctly. As we are working with globally gridded data, we will define `ColumnRingGrid` based on the `global_grid` we already initialized. Then, we will load our inputs. For this, we will choose the January (so the first element) of our climatology files. When using them in `InputSource` be sure to choose the same name and units as used in the definitions of the dynamics before.

snow_grid = ColumnRingGrid(UniformSpacing(N = 1), global_grid, land_sea_mask);
snow_input = InputSource(snow_grid, snow_climatology[:, 1], name = :snow_fall, units = u"m/s");
lst_input = InputSource(snow_grid, lst_climatology[:, 1], name = :air_temperature, units = u"°C")

# As an initial condition, we just cover the whole Earth in deep snow (everywhere the same)!

snow_initializers = (snow_storage = 0.5,)

# Now, we initialize our model and the integrator. As in the first example, we use a `Heun` time stepper

snow_model = SnowModel(snow_grid);
snow_integrator = initialize(snow_model, Heun(Δt = Float32(1.0)), snow_input, lst_input; initializers = snow_initializers)

# ... and we can finally run the model. As before, by wrapping it in an `Oceananigans.Simulation` to output our results
snow_sim = Simulation(snow_integrator; stop_time = 7.0e6, Δt = 3600.0)
Terrarium.initialize!(snow_integrator)
output_dir_snow = mkpath(tempname())
output_file_snow = joinpath(output_dir_snow, "ddsnow-simulation.jld2")
snow_sim.output_writers[:snapshots] = JLD2Writer(
    snow_integrator,
    (snow_storage = snow_integrator.state.snow_storage,);
    filename = output_file_snow,
    overwrite_existing = true,
    including = [:grid],
    schedule = TimeInterval(3600)
)
run!(snow_sim)
@assert isfile(output_file_snow) "Output file does not exist!"
println("Simulation data saved to $(output_file_snow)")

# And now, we plot that data. First we load the JLD2 file.
fts_result = FieldTimeSeries(output_file_snow, "snow_storage")

# Then, we plot it using `CairoMakie`. For this purpose we first convert to a `RingGrids.Field` and then plot it via `heatmap`
tsteps = 1
ring_field = RingGrids.Field(fts_result[tsteps], snow_grid)[:, 1]
heatmap(ring_field)
#
fig = Figure(size = (1200, 660))
ax = Axis(
    fig[1, 1],
    aspect = 2,
    title = "Snow water equivalent [m]",
    titlesize = 20,
    xticks = 0:60:360,
    yticks = -60:30:60,
    xticklabelsize = 10,
    yticklabelsize = 10,
    xtickformat = values -> ["$(round(Int, value))˚E" for value in values],
    ytickformat = values -> ["$(round(Int, value))˚N" for value in values],
)
lond = RingGrids.get_lond(ring_field)
latd = RingGrids.get_latd(ring_field)
n_t = Observable(1)
data = @lift Matrix(RingGrids.Field(fts_result[$n_t], snow_grid)[:, 1])
hm = heatmap!(ax, lond, latd, data, colorrange = (0, 1))
Colorbar(fig[:, end + 1], hm)
frames = 1:length(fts_result)
record(fig, joinpath(@__DIR__, "snow_storage.mp4"), frames, framerate = 12) do i
    n_t[] = i
end

# ![Snow storage animation](snow_storage.mp4)

# And just like that we have implemented our snow simulation. In this version the forcing / input is completly static, so we converge to a static point that corresponds to the January climatology. That's why still see a faily big snow cover in the northern hemisphere. An obvious next step for this model would be now to actually use the full seasonal climatology.
