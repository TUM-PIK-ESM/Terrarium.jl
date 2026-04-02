# Basic concepts

```@meta
CurrentModule = Terrarium
```

Terrarium should not be thought of as a single model but rather a *framework* or "toolkit" for building a wide range of different possible land models from a set of process-based building blocks. This page provides a brief overview of the core abstractions and concepts that form the basis of this framework. 

## Models

A “model” in Terrarium represents a collection of components that fully characterize the dynamics of a land simulation. All models must be defined as `struct`s that subtype [`AbstractModel`](@ref) and typically consist of the following components:
- A `grid` that defines both a vertical and lateral discretization of the spatial domain,
- One or more [Processes](@ref) that define the state variables, parameters, and dynamics of the model,
- An `initializer` that defines a sequence of initialization routines (as well as any associated parameters) for the state variables declared by all of its components.

To see this in action, let's look again at a simple example of setting up a soil model for a single vertical column:

```@example soil_quickstart
using Terrarium

# Set up a SoilModel on a ColumnGrid with 10 exponentially spaced vertical soil layers that will run on the CPU with 32-bit precision
vert = ExponentialSpacing(Δz_min = 0.1, N = 10)
grid = ColumnGrid(CPU(), Float32, vert)
model = SoilModel(grid)
```

Here we first specify the `grid` which defines the spatial discretization (a single rectangular `column`), device architecture, and number format (`Float32`) to be used by all model components. The spacing of the $N = 10$ vertical layers is set to be exponentially increasing from a minimum thickness of 10 cm at the surface; note that alternative options for the vertical discretization include [`UniformSpacing`](@ref) and [`PrescribedSpacing`](@ref). The `SoilModel` is then constructed directly from this `grid` using its default process configurations. Individual processes and parameterizations can also be constructed and passed in as keyword arguments, e.g,

```@example soil_quickstart
hydrology = SoilHydrology(eltype(grid), RichardsEq())
soil = SoilEnergyWaterCarbon(eltype(grid); hydrology)
model = SoilModel(grid; soil)
```

`SoilEnergyWaterCarbon` is a *coupled* process of which `SoilHydrology` is one component. Note that the keyword argument syntax `; hydrology` is a Julia shorthand for `hydrology = hydrology`.

## Processes

Procceses are the building blocks of all Terrarium models. Implementations of `AbstractProcess` represent physical processes characterized by:
- Zero or more state `variable`s that vary spatially across any given `grid`,
- Zero or more *parameters* that are spatially constant and defined somewhere within the process `struct`,
- One or more functions defining equations that compute quantities of interest from both the state variables and parameters defined by the process.

Concrete implementations of `AbstractProcess` are `struct` types that typically consist of zero or more parameters or parameter `struct`s (also sometimes referred to as *parameterizations*). Spatially varying parameters should generally be defined as [`InputVariable`](@ref)s (see following section) or computed by [kernel functions](@ref "Kernel programming") based on spatially invariant parameters stored in parameterization `struct`s.

## State variables
State variables are symbolically defined for each process in their respective implementations of the `variables` method by returning instances of [`AbstractVariable`](@ref) which define (at minimum) its name, dimensionality, and physical units. State variables may be one of three types:
```@docs; canonical = false
PrognosticVariable
AuxiliaryVariable
InputVariable
```

These variables are typically constructed using one of the convenience methods [`prognostic`](@ref), [`auxiliary`](@ref), or [`input`](@ref).

!!! info "Variables vs. `Field`s"
    Within the context of defining and handling model/process state, Terrarium distinguishes between *variables* and *Fields*. Variables are symbolic objects subtyping [`AbstractVariable`](@ref) that contain metadata about the state variable defined by the process or model type. These variables are realized as [`Fields`](@ref) defined over the model grid that contain the actual numerical data of the corresponding variable. 

As a simple example, suppose we are implementing a new process `MyProcess` and we want to define the necessary state variables. We do this by defining a new dispatch of the `variables` method:
```julia
variables(::MyProcess) = (
    prognostic(:progvar, XYZ()),
    auxiliary(:auxvar, XYZ()),
    auxiliary(:bc, XY()),
    input(:input, XY())
)
```
This will result in a total of five state variables being allocated upon initialization: two auxiliary variables named `auxvar` and `bc` as well as one prognostic variable named `progvar` along with its corresponding tendency variable which is created automatically. The second argument to the variable metadata constructors `prognostic` and `auxiliary` is a subtype of `VarDims` which specifies on which spatial dimensions the state variable should be defined. `XYZ()` corresponds to a 3D `Field` which varies both laterally and with depth.

Currently, Terrarium only supports a single 3D grid representing variables defined in the soil domain, though this may change in the future in order to accommodate multi-layer snow and canopy processes. `XY()` corresponds to a 2D field which is discretized along the lateral dimension only. Note that Terrarium also currently supports only 1D (vertical) dynamics so all grid cells on the X and Y axes will be assumed independent. This is equivalent to what is typically called a *single column* model, or *column-based parameterization* in atmosphere and ocean modeling. However, building on Oceananigans means that we have a clear path to relax this assumption in the future!
