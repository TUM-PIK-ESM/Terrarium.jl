# [Vegetation](@id vegetation_docs)

```@meta
CurrentModule = Terrarium
```

```@setup vegetation
using Terrarium
using InteractiveUtils
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

Vegetation processes simulate the biological and biogeochemical dynamics of terrestrial plants, including photosynthesis, phenology, respiration, and carbon cycling. The vegetation module in Terrarium currently consists of the following processes:

- **Photosynthesis** — light-dependent gross primary production (CO₂ uptake)
- **Stomatal conductance** — controls leaf gas exchange in response to light, humidity, and water stress
- **Autotrophic respiration** — plant maintenance and growth respiration
- **Phenology** — seasonal variation in leaf area (LAI) and fractional vegetation coverage
- **Carbon dynamics** — prognostic evolution of vegetation carbon pools (leaves, stems, roots, etc.)
- **Plant-available water** (PAW) — soil moisture stress factor affecting transpiration and photosynthesis
- **Root distribution** (optional) — vertical profile of root density for soil-plant water coupling
- **Vegetation dynamics** (optional) — prognostic plant population density or fractional coverage evolution

These individual processes are composed together in implementations [`AbstractVegetation`](@ref) which provide complete representations of vegetation carbon and population dynamics with which other model components can be coupled.

```@docs; canonical = false
AbstractVegetation
```

```@example vegetation
subtypes(Terrarium.AbstractVegetation)
```

## Vegetation carbon

Terrarium currently provides a single implementation of [`AbstractVegetation`](@ref), [`VegetationCarbon`](@ref), which couples together all of the above processes related to the carbon cycle for natural vegetation. `VegetationCarbon` also provides a coupling interface for interacting with [`AbstractAtmosphere`](@ref) and [`AbstractSoil`](@ref) components.

```@docs; canonical = false
VegetationCarbon
```

```@example vegetation
variables(VegetationCarbon(Float32))
```

### Process interface

```@docs; canonical = false
compute_auxiliary!(
        state, grid,
        veg::VegetationCarbon,
        constants::PhysicalConstants,
        atmos::AbstractAtmosphere,
        soil::Optional{AbstractSoil} = nothing,
        args...
    )

compute_tendencies!(state, grid, veg::VegetationCarbon, args...)
```

## Component processes

Detailed implementations, parameterization options, and kernel functions for individual vegetation processes are documented on dedicated pages:

- [Photosynthesis](@ref) — photosynthesis and CO₂ uptake
- [Stomatal conductance](@ref) — stomatal regulation and gas exchange
- [Autotrophic respiration](@ref) — plant respiration rates
- [Phenology](@ref) — LAI dynamics and seasonal cycles
- [Plant available water](@ref) — soil moisture stress parameterizations
- [Vegetation carbon dynamics](@ref) — carbon pool evolution and allocation
- [Root distribution](@ref) — vertical root profiles and soil-plant coupling
