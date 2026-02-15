# Baisc concepts

Terrarium should not be thought of as a single model but rather a *framework* or "toolkit" for building a wide range of different possible land models from a set of process-based building blocks. This page provides a brief overview of the core abstractions and concepts that form the basis of this framework. 

## Models

A “model” in Terrarium represents a collection of components that fully characterize the dynamics of a land simulation. All models must be defined as `struct`s that subtype [`AbstractModel`](@ref) and typically need consist of the following components:
- A `grid` that defines both a vertical and lateral discretization of the spatial domain,
- One or more [Processes](@ref) that define the state variables, parameters, and dynamics of the model,
- An `initializer` that defines a sequence of initialization routines (as well as any assocaited parameters) for the state variables declared by all of its components.

To see this in action, let's look again at a simple example of setting up a soil model for a single vertical column:

```@example soil_quickstart
using Terrarium

# Set up a SoilModel on a ColumnGrid with 10 exponentially spaced vertical soil layers that will run on the CPU with 32-bit precision
vert = ExponentialSpacing(Δz_min = 0.1, N = 10)
grid = ColumnGrid(CPU(), Float32, vert)
model = SoilModel(grid)
```

Here we first specify the `grid` which defines the spatial discretization (a single rectangular `column`), device architecture, and number format (`Float32`) to be used by all model components. The spacing of the $N = 10$ vertical layers is set to be exponentially increasing from a minimum thickness of 10 cm at the surface; note that alternative options for the vertical discretization include [`UniformSpacing`](@ref) and [`PrescribedSpacing`](@ref).

The `SoilModel` is then constructed directly from this `grid`. If the 

## Processes

Implementations of `AbstractProcess` represent physical processes characterized by:
- Zero or more state `variable`s that vary spatially across any given `grid`,
- Zero or more *parameters* that are spatially constant and defined somewhere within the process `struct`,
- One or more functions defining equations that compute quantities of interest from both the state variables and parameters defined by the process.

Concrete implementations of `AbstractProcess` are `struct` types that typically consist of zero or more parameters or parameter `struct`s (also sometimes referred to as *parameterizations*).

### Parameterizations

## State variables

## Initialization

## Timestepping
