# Initialization

```@meta
CurrentModule = Terrarium
```

## Overview

Model and process types in Terrarium are *stateless* and *immutable*; i.e. they only specify parameters and model configuration. To allocate state variables for a model or process, we need to use the [`initialize`](@ref) method:

```@docs; canonical = false
initialize(model::AbstractModel{NF, Grid}) where {NF, Grid}
initialize(process::AbstractProcess{NF}, grid::AbstractLandGrid{NF}) where {NF}
```

which will create and return a [`StateVariables`](@ref) structure containing all of the initialized [`Field`](@ref)s corresponding to state variables defined by the model/process.

```@docs; canonical = false
StateVariables
```

By default, `Field`s will be initialized with zeros. Some state variable types, such as `input` variables, allow for the specification of alternative default values.

Of course, very few models can do anything useful or interesting starting from (literally) zero. Terrarium provides three complementary mechanisms through which initialization routines for model/process state variables can be defined:

- Direct initialization via the `initializers` keyword argument of [`initialize`](@ref)
- [`AbstractInitializer`](@ref) types specified in models encapsulate a sequence of initialization routines
- Model and process specific dispatches of `initialize!`

As a general rule, these initializers are invoked in the order that they are listed above. Implementations of `initialize!` should generally assume that the prognostic state has already been initialized by corresponding `Field` or model initializers.

## Direct initialization of `Field`s

The keyword argument must be a `NamedTuple` where the keys correspond to the name of the state variable and the values are either scalars, arrays matching the size of the model `grid`, or functions of the form `f(coords...)` where `coords` are the non-[`Flat`](@extref Oceananigans.Grids.Flat) dimensions of `grid`. For column-based grids, this is generally `f(x,z)` with `x` corresponding to a column index.

## Model initializers

These `Initializer` types can be supplied to subtypes of [`AbstractModel`](@ref) during construction. Models can/should typically define corresponding [`AbstractInitializer`](@ref) types that represent common initialization strategies appropriate for the processes included in that model; e.g. the [`SoilModel`](@ref) defines [`SoilInitializer`](@ref) with process-specific initialization types like [`QuasiThermalSteadyState`](@ref) and [`SaturationWaterTable`](@ref).

## Built-in initialization routines

Model/process specific dispatches of [`initialize!`](@ref) can be defined for process initialization logic that should be run regardless of the choice of prognostic state variable initialization.
