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
begin
	using Terrarium 
	
	@kwdef struct ExpModel{NF, Grid, I, BC} <: Terrarium.AbstractModel{NF, Grid}
	    grid::Grid
	    initializer::I = DefaultInitializer()
	    boundary_conditions::BC = DefaultBoundaryConditions()
	end
	
	ExpModel(grid::Grid, initializer::I, boundary_conditions::BC) where {Grid, I, BC} = ExpModel{eltype(grid), Grid, I, BC}(grid, initializer, boundary_conditions)
end

# ╔═╡ 5630efd5-2482-463d-913f-9addb120beec
md"""
# A super basic example model with Terrarium 

In this example we will set up an embarrassingly simple example to demonstrate Terrarium's model interface. Our model should just exhibit 1-dimensional exponential dynamics with an offset as follows

```math 
\frac{du}{dt} = u + c
```

with some state variable ``u``. For the sake of this demonstration we will treat the offset ``c`` as an auxiliary/diagnostic variable even though it is constant in time.

We begin by defining our model `struct` that subtypes `Terrarium.AbstractModel`: 
"""

# ╔═╡ 4922e264-c80d-4a5b-8891-a7c8a3fdbfe7
md""" 
A `Terrarium.AbstractModel` typically constists of the fields 
 * `grid` which defines the discretization of the spatial domain
 * `initializer` which is responsible for initializing state variables
 * `boundary_conditions` that defines the behaviour fo the state variabels at the spatial boundaries of the `grid` 
 * further fields that define processes, dynamics and submodels 

When we follow the advised naming notations of `grid`, `initializer` and `boundary_conditions` we inherit default methods from `Terrarium.AbstractModel` that we will suffice for our basic example. For more complex models we might need to implement custom `initialize!(state, ::Model, ::Initializer)` to initialize model states. 

## What's the `grid`? 

The `grid` defines the spatial discretization. Our implementation of the grid uses Oceananigans.jl (and SpeedyWeather.jl/RingGrids.jl) grids to ease computations (and coupling). We have two main grid types available: 

* `ColumnGrid` is a set of laterally independent vertical columns with dimensions ``(x, y, z)`` where ``x`` is the column dimension, ``y=1`` is constant, and ``z`` is the vertical axis, 
* `ColumnRingGrid` represents a global (spherical) grid of independent, vertical columns where the spatial discretization in the horizontal direction is defined by a `RingGrids.AbstractGrid`. 

In both cases we need to specificy the vertical discretizataion via an `UniformSpacing`, `ExponentialSpacing` or `PrescribedSpacing`

## Initializer and Boundary Conditions 

For our basic example here the defaults will suffice, and we won't have to define custom ones. 

## What's our `grid`? 

For our example we just want a single column with a single level, we can define it like so: 
"""

# ╔═╡ 78f268ef-5385-4c63-bc35-2c973de69da5
grid = ColumnGrid(CPU(), Float64, UniformSpacing(N=1))

# ╔═╡ 575d920c-b12e-493f-95a7-5c962c3591fd
md"""
## Defining the model behaviour 

Now, we want to define our intended model behaviour. For this, we need to define the following routines: 

* `variables(::Model)` returns a tuple of variable metadata declaring the state variables. Variables must be one of two types: `prognostic` or `auxiliary` (sometimes referred to as “diagnostic”). Prognostic variables fully characterize the state of the system at any given timestep and are updated according to their tendencies (i.e. ``u`` in our example). Tendencies are automatically allocated for each prognostic variable declared by the model. In this example we will treat the offset ``c`` as an auxiliary variable (even though we could just include it as a constant in the tendency computations)
* `compute_auxiliary!(state, ::Model)` computes the values of all auxiliary variables (if necessary) assuming that the prognostic variables of the system in state are available for the current timestep.
* `compute_tendencies!(state, ::Model)` computes the tendencies based on the current values of the prognostic and auxiliary variables stored in state.

So, let's define those: 
"""

# ╔═╡ 82e45724-ba16-4806-9470-5cb4c43ea734
Terrarium.variables(::ExpModel) = (Terrarium.prognostic(:u, Terrarium.XY()), 
                              	   Terrarium.auxiliary(:c, Terrarium.XY()))

# ╔═╡ 32373599-768f-4809-acdd-4704acc3f30b
md"""
Here, we defined our two variables with their name as a `Symbol` and whether they are 2D variables (`XY`) or 3D variables (`XYZ`). 

## Our state variables

Next, we define the actual compute function. For this we need to know a little about how the variables we just defined are handled in our `StateVariables`. The `StateVariables` hold all prognostic and auxiliary variables, their tendencies and closures and additional inputs and forcings in seperate `NamedTuples`. We have shortcuts defined so, that e.g. in our example both `state.prognostic.u` and `state.u` will work. With that let's define the following:
"""

# ╔═╡ d55aaf4c-3033-45ba-9d64-8fa8ae4b671c
function Terrarium.compute_auxiliary!(state, model::ExpModel) 
    state.auxiliary.c .= 0.1
end 

# ╔═╡ 5c8be7e4-f150-492b-a75d-96887a11f6da
# du/dt = u + c = u + 0.1
function Terrarium.compute_tendencies!(state, model::ExpModel) 
    state.tendencies.u .= state.prognostic.u + state.auxiliary.c
end

# ╔═╡ 8d331856-6e9b-41d4-b1be-a84d5fedac8d
md"""
These example compute functions are really the simplest possible, for more complex operations, we would need to define them via `KernelAbstractions` kernels. In this notebook we will not go into further details on that. 

However, now we have everything our model needs and we can use it! 

## Running our model 

First, we need to define our initial conditions via a `FieldInitializer`.

`FieldInitializer` can take in functions `(x,z)->val`, arrays or values. It uses `Oceananigans.set!(field, x)`, and allows all input arguments for `x` that `set!` allows: 
"""

# ╔═╡ fad517e2-9d3a-43cd-8e4d-0b72ca55b8c8
initializer = FieldInitializers(u = 0., c = 0.1)

# ╔═╡ 4b483c23-9e15-4d03-b275-8f530854669e
md"""
Then, we can actually construct our model 
"""

# ╔═╡ 2a4234c5-f529-4166-94c3-0556565348ea
model = ExpModel(; grid, initializer)

# ╔═╡ 4c36fdc0-5120-46b9-86ca-e875e23a6c1d
md"""
Then we initialize our model, i.e. we run all pre-computation, and initialize the `StateVariables`. Here, we can also define our timestepper
"""

# ╔═╡ 7e38132b-d406-4863-b88f-90efe2a1bfa2
modelstate = initialize(model; timestepper=Heun)

# ╔═╡ 21e20c28-dfe1-4a0a-992f-c3499fbe4be8
md"""
Then, we can finally run our model using `run!` that allows either for a number of `steps` or a `Dates.Period` to be specified: 
"""

# ╔═╡ de3d4210-c39f-11f0-3d50-3f95a2361e2a
run!(modelstate, period=Hour(1))

# ╔═╡ cce4d4d3-0fa4-4376-bcb6-c52603bc17d6
modelstate.state.u

# ╔═╡ 7fa2dfbf-7077-4162-bcc1-ba2bd12b093c
md"""
In one hour our state already grew to `7.46953e54`, if that's not exponential growth ;) 

Well, and that's already it. We defined a simple exponential model within our Terrarium `AbstractModel` interface!
"""

# ╔═╡ Cell order:
# ╟─5630efd5-2482-463d-913f-9addb120beec
# ╠═808d5d89-c1d2-4f6a-bd55-4b3a8444c90f
# ╠═94d82d31-42ec-41de-91e9-b5585b3a72d4
# ╟─4922e264-c80d-4a5b-8891-a7c8a3fdbfe7
# ╠═78f268ef-5385-4c63-bc35-2c973de69da5
# ╟─575d920c-b12e-493f-95a7-5c962c3591fd
# ╠═82e45724-ba16-4806-9470-5cb4c43ea734
# ╟─32373599-768f-4809-acdd-4704acc3f30b
# ╠═d55aaf4c-3033-45ba-9d64-8fa8ae4b671c
# ╠═5c8be7e4-f150-492b-a75d-96887a11f6da
# ╟─8d331856-6e9b-41d4-b1be-a84d5fedac8d
# ╠═fad517e2-9d3a-43cd-8e4d-0b72ca55b8c8
# ╟─4b483c23-9e15-4d03-b275-8f530854669e
# ╠═2a4234c5-f529-4166-94c3-0556565348ea
# ╟─4c36fdc0-5120-46b9-86ca-e875e23a6c1d
# ╠═7e38132b-d406-4863-b88f-90efe2a1bfa2
# ╟─21e20c28-dfe1-4a0a-992f-c3499fbe4be8
# ╠═de3d4210-c39f-11f0-3d50-3f95a2361e2a
# ╠═cce4d4d3-0fa4-4376-bcb6-c52603bc17d6
# ╟─7fa2dfbf-7077-4162-bcc1-ba2bd12b093c
