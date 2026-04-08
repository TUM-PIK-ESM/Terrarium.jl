# Land model

```@meta
CurrentModule = Terrarium
```

```@setup landmodel
using Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

[`LandModel`](@ref) is a fully-coupled land surface and terrestrial ecosystem model that simulates the integrated dynamics of soil, surface hydrology, surface energy balance, vegetation, and atmospheric forcing. Conceptually, the model represents a single, homogeneous land surface tile within each `grid` cell with a predetermined set of processes.

```@example landmodel
arch = CPU()
grid = ColumnGrid(arch, Float32, ExponentialSpacing(N = 10)) # 10 soil layers
model = LandModel(grid) # Default configuration
integrator = initialize(model, ForwardEuler(eltype(grid)))
```

```@docs; canonical = false
LandModel
```

```@example landmodel
variable(model)
```

## Components

[`LandModel`](@ref) integrates five major sub-processes: **atmosphere**, **soil**, **surface energy balance**, **surface hydrology**, and **vegetation**. Each component can be configured separately when constructing a `LandModel`; see the linked process pages for available implementations and parameterizations.

| Field | Type | Scope | Process page |
|-------|------|-------|---------------|
| `atmosphere` | [`AbstractAtmosphere`](@ref) | Meteorological input variables | [Atmosphere](@ref) |
| `soil` | [`AbstractSoil`](@ref) | Energy, water, carbon in soil | [Soil processes](@ref) |
| `surface_energy_balance` | [`AbstractSurfaceEnergyBalance`](@ref) | Radiative and turbulent energy fluxes | [Surface energy balance](@ref) |
| `surface_hydrology` | [`AbstractSurfaceHydrology`](@ref) | Infiltration, evapotranspiration, interception | [Surface hydrology](@ref) |
| `vegetation` | `Optional{`[`AbstractVegetation`](@ref)`}` | Coupled vegetation carbon processes | [Vegetation](@ref) |

Each component is summarized briefly below. See the linked pages for further details about each component process.

### Atmosphere

The `atmosphere` component provides (possibly time-varying) meteorological inputs. The default implementation is [`PrescribedAtmosphere`](@ref), which reads air temperature, humidity, wind, radiation, and precipitation from [`InputVariable`](@ref)s and provides them as boundary conditions to the surface energy balance, hydrology, and vegetation components. See [Atmospheric inputs](@ref) for further details on the atmosphere interface.

### Soil

The `soil` component represents the solid land surface extending from the topmost soil layer down to an arbitrary depth determined by the `grid`. The default configuration of `LandModel` uses [`SoilEnergyWaterCarbon`](@ref) which represents coupled energy, water, and carbon transport within the soil column. See [Soil processes](@ref) for detailed descriptions of energy, hydrology, and biogeochemistry implementations.


## Initializers

!!! todo "Initializers and boundary conditions"
    `LandModel` does not yet have its own dedicated boundary conditions and initializers, but it will soon! Stay tuned!
