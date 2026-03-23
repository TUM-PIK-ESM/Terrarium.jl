# Soil model

```@meta
CurrentModule = Terrarium
```

```@setup default
using Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

[`SoilModel`](@ref) is a general 1D column model that couples soil energy, water, and carbon
dynamics. It is the primary model type for simulating heat conduction, freeze-thaw processes,
and variably saturated water flow in the subsurface. All prognostic state variables are
discretized in space and time via user-specified `grid` and time stepper.

```@example default
arch = CPU()
grid = ColumnGrid(arch, Float32, ExponentialSpacing(N = 10)) # 10 soil layers
model = SoilModel(grid)
integrator = initialize(model, ForwardEuler(eltype(grid)))
```

## Abstract types

```@docs; canonical = false
AbstractSoilModel
```

```@docs; canonical = false
AbstractSoil
```

## Concrete types

```@docs; canonical = false
SoilModel
```

### Components

The default soil process bundle is [`SoilEnergyWaterCarbon`](@ref), which couples four
sub-processes. All four can be replaced independently when constructing a `SoilModel`; see
each linked process page for available implementations and parameterizations.

```@docs; canonical = false
SoilEnergyWaterCarbon
```

| Field | Type | Process page |
|-------|------|--------------|
| `strat` | [`AbstractStratigraphy`](@ref) | [Soil stratigraphy](@ref) |
| `energy` | [`AbstractSoilEnergyBalance`](@ref) | [Soil energy balance](@ref) |
| `hydrology` | [`AbstractSoilHydrology`](@ref) | [Soil hydrology](@ref) |
| `biogeochem` | [`AbstractSoilBiogeochemistry`](@ref) | — |

Each component is summarized briefly below. Follow the linked process pages for full
theoretical background, available concrete types, state variables, and method signatures.

#### Stratigraphy

The `strat` field parameterizes the vertical distribution of soil material properties (texture,
porosity, organic content). It provides kernel functions used by the energy and hydrology
sub-processes to look up spatially varying material properties at each grid cell. By default
[`HomogeneousStratigraphy`](@ref) is used, which assumes a single uniform material throughout
the profile. See [Soil stratigraphy](@ref) for details.

#### Energy balance

The `energy` field represents heat conduction in the soil column, including the latent heat of
freeze-thaw phase change. The default implementation is [`SoilEnergyBalance`](@ref), which
evolves the volumetric internal energy $U$ (J m⁻³) as the prognostic variable and derives
temperature via the [`SoilEnergyTemperatureClosure`](@ref). See [Soil energy balance](@ref) for
details.

#### Hydrology

The `hydrology` field governs the vertical movement of water in the soil column. The default
implementation is [`SoilHydrology`](@ref), which evolves total saturation (liquid + ice) as
the prognostic variable. Vertical flow is disabled by default ([`NoFlow`](@ref)); enabling
the Richards equation requires selecting [`RichardsEq`](@ref) as the `vertflow` operator. See
[Soil hydrology](@ref) for details.

#### Biogeochemistry

The `biogeochem` field provides the spatial distribution of soil organic carbon and associated
biogeochemical fluxes. The default implementation [`ConstantSoilCarbonDensity`](@ref) prescribes
a spatially homogeneous organic carbon density and does not evolve any prognostic variables.

## Initialization

Terrarium provides a hierarchy of initializers for `SoilModel`. The top-level initializer is
[`SoilInitializer`](@ref), which composes separate sub-initializers for energy, hydrology, and
biogeochemistry. It is the recommended initializer for most use cases:

```@docs; canonical = false
SoilInitializer
```

Passing a `SoilInitializer` to `SoilModel` overrides the `DefaultInitializer`:

```@example default
initializer = SoilInitializer(Float32;
    energy    = QuasiThermalSteadyState(Float32; T₀ = 2.0f0),
    hydrology = SaturationWaterTable(Float32; water_table_depth = 3.0f0),
)
model = SoilModel(grid; initializer)
```

### Energy initializers

```@docs; canonical = false
ConstantInitialSoilTemperature
```

```@docs; canonical = false
QuasiThermalSteadyState
```

```@docs; canonical = false
PiecewiseLinearInitialSoilTemperature
```

### Hydrology initializers

```@docs; canonical = false
SaturationWaterTable
```

### Fallback

```@docs; canonical = false
DefaultInitializer
```

## Boundary conditions

Terrarium provides a set of named boundary condition aliases for the most common `SoilModel`
configurations. These return `NamedTuple`s in the format expected by the `boundary_conditions`
keyword argument of [`initialize`](@ref).

### Energy boundary conditions

```@docs; canonical = false
GroundHeatFlux
```

```@docs; canonical = false
GeothermalHeatFlux
```

```@docs; canonical = false
PrescribedSurfaceTemperature
```

```@docs; canonical = false
PrescribedBottomTemperature
```

### Hydrology boundary conditions

```@docs; canonical = false
InfiltrationFlux
```

```@docs; canonical = false
ImpermeableBoundary
```

```@docs; canonical = false
FreeDrainage
```
