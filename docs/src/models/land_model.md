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

[`LandModel`](@ref) is a fully-coupled land surface and terrestrial ecosystem model that simulates the integrated dynamics of soil, surface hydrology, surface energy balance, vegetation, and atmospheric forcing. Conceptually, [`LandModel`](@ref) represents a single, homogeneous land surface tile within each `grid` cell with a predetermined set of processes.

```@example landmodel
arch = CPU()
grid = ColumnGrid(arch, Float32, ExponentialSpacing(N = 10)) # 10 soil layers
model = LandModel(grid) # Default configuration
integrator = initialize(model, ForwardEuler(eltype(grid)))
```

```@docs; canonical = false
LandModel
```

## Components

[`LandModel`](@ref) integrates five major sub-processes: **atmosphere**, **soil**, **surface energy balance**, **surface hydrology**, and **vegetation**. Each component can be configured separately when constructing a `LandModel`; see the linked process pages for available implementations and parameterizations.

| Field | Type | Scope | Process page |
|-------|------|-------|---------------|
| `atmosphere` | [`AbstractAtmosphere`](@ref) | Meteorological and orbital input variables | [Atmosphere](@ref) |
| `soil` | [`AbstractSoil`](@ref) | Energy, water, carbon in soil | [Soil processes](@ref) |
| `surface_energy_balance` | [`AbstractSurfaceEnergyBalance`](@ref) | Radiative and turbulent energy fluxes | [Surface energy balance](@ref) |
| `surface_hydrology` | [`AbstractSurfaceHydrology`](@ref) | Infiltration, evapotranspiration, interception | [Surface hydrology](@ref) |
| `vegetation` | `Optional{`[`AbstractVegetation`](@ref)`}` | Plant phenology, photosynthesis, carbon | [Vegetation](@ref) |

Each component is summarized briefly below. Follow the linked process pages for full theoretical background, available concrete types, state variables, and method signatures.

### Atmosphere

The `atmosphere` component provides (possibly time-varying) meteorological inputs. The default implementation is [`PrescribedAtmosphere`](@ref), which reads air temperature, humidity, wind, radiation, and precipitation from [`InputVariable`](@ref)s and provides them as boundary conditions to the surface energy balance, hydrology, and vegetation components. See [Atmospheric inputs](@ref) for further details on the atmosphere interface.

### Soil

The `soil` component represents the solid land surface extending from the topmost soil layer down to an arbitrary depth determined by the `grid`. The default configuration of `LandModel` uses [`SoilEnergyWaterCarbon`](@ref) which represents coupled energy, water, and carbon transport within the soil column. See [Soil processes](@ref) for detailed descriptions of energy, hydrology, and biogeochemistry implementations.

### Surface Energy Balance

The `surface_energy_balance` component calculates the energy budget at the land surface, balancing net radiation with sensible, latent, and ground heat fluxes. The default implementation is [`SurfaceEnergyBalance`](@ref), which couples radiative and turbulent flux closures and prognostic or diagnostic skin temperature. See [Surface energy balance](@ref) for available radiative and turbulent parameterizations.

### Surface Hydrology

The `surface_hydrology` component partitions mass and energy fluxes at the surface via canopy interception, evapotranspiration, and infiltration/runoff. The default is [`SurfaceHydrology`](@ref), which couples these processes and uses soil moisture availability to constrain evapotranspiration. See [Surface hydrology](@ref) for available ET schemes and interception/runoff parameterizations.

### Vegetation

The `vegetation` component (optional; default: [`VegetationCarbon`](@ref)) represents the terrestrial biosphere, including photosynthesis, phenology, carbon cycling, and plant-available water. Vegetation modulates the surface by absorbing radiation, controlling evapotranspiration via stomata, and extracting soil water through roots. Set `vegetation = nothing` for bare-ground simulations. See [Vegetation](@ref) for an overview of the relevant processes and parameterizations.

## Initializers

!!! todo "Initializers and boundary conditions"
    `LandModel` does not yet have its own dedicated boundary conditions and initializers, but it will soon! Stay tuned!
