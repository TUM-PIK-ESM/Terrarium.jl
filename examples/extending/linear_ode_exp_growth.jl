# # Getting started with a simple exponential growth model
#
# In this example, we will set up an embarrassingly simple example to demonstrate Terrarium's model interface. Our model will have 1-dimensional exponential dynamics with a constant offset
#
# ```math
# \frac{du}{dt} = \alpha u + c + F(t)
# ```
#
# for an arbitrary prognostic variable ``u``. For the sake of this demonstration we will treat the offset ``c`` as an auxiliary/diagnostic variable even though it is constant in time. ``F(t)`` is an external forcing that we apply.

using Terrarium

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
# In both cases we need to specify the vertical discretization via an `UniformSpacing`, `ExponentialSpacing` or `PrescribedSpacing`.
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

# ## Defining the model behavior
#
# Now, we want to define our intended model behavior. For this, we need to define the following methods:
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

using Random

Random.seed!(1234) # set random seed

t_F = 0:1:300;
F = FieldTimeSeries(grid, XY(), t_F);
F.data .= randn(size(F));
input = InputSource(grid, F, name = :F)

# Here we constructed a 2D (`XY()`) time series on our `grid` at times `t_F` with random normal distributed data and defined our `InputSource` for our model based on it.

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

using Oceananigans: JLD2Writer, TimeInterval
using Oceananigans.Units: seconds
using JLD2

Terrarium.initialize!(integrator)
output_dir = mkpath(tempname())
output_file = joinpath(output_dir, "simulation.jld2")
sim.output_writers[:snapshots] = JLD2Writer(
    integrator,
    (u = integrator.state.u,);
    filename = output_file,
    overwrite_existing = true,
    including = [:grid], # include the grid with the output
    schedule = TimeInterval(10seconds)
)
run!(sim)
@assert isfile(output_file) "Output file does not exist!"
println("Simulation data saved to $(output_file)")

# Then load the output data and plot the results:

using CairoMakie

fts = FieldTimeSeries(output_file, "u")

plot(1:length(fts), [fts[i][1, 1, 1] for i in 1:length(fts)])
