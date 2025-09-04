# Software architecture

!!! warning "🚧🚧 Under construction 🚧🚧"
    The architecture of Terrarium.jl is still in an early prototype stage and is the subject of ongoing discussion. This page is a non-exhaustive summary of the current working concept and will be updated accordingly if and when the architecture changes.

## Numerical implementation
Terrarium is based on the numerics and finite-volume method operators provided by Oceananigans.jl. All state variables are represented as Oceananigans Fields with computations automatically distributed using KernelAbstractions.jl. Spatial grids in Terrarium are similarly based on the corresponding Oceananigans grids, most prominently the RectilinearGrid which represents a rectangular volume divided orthogonally into smaller control volumes along the X, Y, and Z dimensions.

## The AbstractModel interface
A “model” in Terrarium represents a collection of parameters and process types that characterize a simulation. All models must be defined as `struct`s that subtype `AbstractModel` and will typically consist of the following fields/properties:
- A grid that defines the discretization of the spatial domain
- An initializer responsible for initializing state variables
- One or more processes and/or sub-models that define the state variables and the dynamics of the system

The distinction between what constitutes a “model” versus a “process” should be based on whether or not there is a use case for running simulations of that model/process independently from other processes. Models should be built to be independently runnable standalone in a `Simulation`, while process types simply provide the concrete implementations of the physical processes to be included in a model. Since Terrarium models are composable, it should generally be trivial to convert a process into a full standalone model later (if necessary) without affecting the other model components that depend on it.

All `AbstractModel` and `AbstractProcess` types must additional provide dispatches for the following methods:
- `variables(::Model)` returns a tuple of variable metadata declaring the state variables. Variables must be one of two types: prognostic or auxiliary (sometimes referred to as “diagnostic”). Prognostic variables fully characterize the state of the system at any given timestep and are updated according to their tendencies (i.e. $G$ in the aforementioned equation). Tendencies are automatically allocated for each prognostic variable declared by the model.
- `compute_auxiliary!(state, ::Model)` computes the values of all auxiliary variables (if necessary) assuming that the prognostic variables of the system in `state` are available for the current timestep.
- `compute_tendencies!(state, ::Model)` computes the tendencies based on the current values of the prognostic and auxiliary variables stored in `state`.

In addition, `AbstractModel`  implementations must also provide dispatches for the following methods:
- `initialize!(state, ::Model, ::Initializer)` computes any necessary initialization of the model state based on the user-supplied configuration and parameter settings. The additional initializer argument is extracted from the model and allows for alternative dispatches based on various initialization schemes for each model.
- `timestep!(state, ::Model, ::TimeStepper, dt)` updates the prognostic state variables according to the given timestepping scheme. This method can and should be implemented generically for any timestepping scheme but allows for model-specific overrides where necessary.
- `get_grid(::Model)` returns the spatial grid associated with the model. A default implementation is provided which assumes that the model defines a field named `grid`.
- `get_initializer(::Model)` returns the initializer for the model. A default implementation is provided which assumes that the model defines a field named `initializer`.

## State variables
State variables are realized as Oceananigans `Field`s defined over the model grid. Each variable returned by the `variables` method will be allocated a corresponding `Field` with boundary conditions and initial values specified by one or more `AbstractInitializer`s defined on the model. This allocation occurs when creating a `Simulation` by calling  `initialize(model)`. The resulting `Simulation` object will have a property `state` of type `StateVariables` which holds all of the `Field`s for each state variable defined by the model.

As an example, suppose we are implementing a new model `MyModel` and we want to define the necessary state variables. We do this by defining a new dispatch of the `variables` method:
```julia
variables(::MyModel) = (
	prognostic(:progvar, XYZ()),
	auxiliary(:auxvar, XYZ()),
	auxiliary(:bc, XY()),
)
```
This will result in a total of four state variables being allocated when `initialize(model)` is called: two auxiliary variables named `auxvar` and `bc` as well as one prognostic variable named `progvar` along with its tendency `progvar_tendency` which is created automatically. The second argument to the variable metadata constructors `prognostic` and `auxiliary` is a subtype of `VarDims` which specifies on which spatial dimensions the state variable should be defined. `XYZ()` corresponds to a 3D `Field` which is discretized both laterally (along spatial grid cells) and along the vertical axis. Currently, Terrarium only supports a single 3D grid representing variables defined in the soil domain, though this may change in the near future in order to accommodate multi-layer snow and canopy processes. `XY()` corresponds to a 2D field which is discretized along the lateral dimension only. Note that Terrarium also currently supports only 1D (vertical) dynamics so all grid cells on the X and Y axes will be assumed independent. This is equivalent to what is typically called a *single column* model, or *column-based parameterization* in atmosphere and ocean modeling. However, building on Oceananigans means that we have a clear path to relax this assumption in the future!
