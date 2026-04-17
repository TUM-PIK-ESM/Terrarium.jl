# Plant available water

```@meta
CurrentModule = Terrarium
```

```@setup vegpaw
using Terrarium
using InteractiveUtils
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

The plant available water (PAW) represents the fraction of soil water that is available for uptake by plant roots. This quantity is critical for constraining photosynthesis and transpiration fluxes.

```@docs; canonical = false
AbstractPlantAvailableWater
```

```@example vegpaw
subtypes(Terrarium.AbstractPlantAvailableWater)
```

### Soil moisture limiting factor

The impact of soil water availability on photosynthesis is aggregated across the rooting depth using the [root distribution](@ref "Root distribution"):
```math
\begin{equation}
\beta_s = \int_0^{z_{\text{max}}} W(z) \cdot r(z) \, dz
\end{equation}
```

where $r(z)$ is the normalized root fraction at depth $z$, and $z_{\text{max}}$ is the rooting depth. This soil moisture limiting factor $\beta_s \in \left[0, 1\right]$ directly modulates the rate of [autotrophic respiration](@ref "Autotrophic respiration") and thus also the net primary plant productivity of vegetation.

## Implementations

### Field capacity limited PAW

A common assumption is to consider PAW to be the fraction of water between the wilting point (where plants can no longer extract water) and field capacity (approximately optimal water availability). The fraction of available water for a given soil layer is then defined as,
```math
\begin{equation}
W = \min\left(\frac{\theta_{\text{liq}} - \theta_{\text{wp}}}{\theta_{\text{fc}} - \theta_{\text{wp}}}, 1\right)
\end{equation}
```
where $\theta_{\text{liq}}$ is the current volumetric water content, $\theta_{\text{wp}}$ is the wilting point, and $\theta_{\text{fc}}$ is the field capacity matric potential. Both wilting point and field capacity are derived from soil texture properties and the water retention characteristics of the soil (see [the documentation on soil hydraulic properties](@ref "Hydraulic properties")).

```@docs; canonical = false
FieldCapacityLimitedPAW
```

```@example vegpaw
variables(FieldCapacityLimitedPAW(Float32))
```

## Process interface

```@docs; canonical = false
compute_auxiliary!(state, grid, paw::FieldCapacityLimitedPAW, soil::AbstractSoil, args...)
```

## Methods

```@docs; canonical = false
soil_moisture_limiting_factor
```

## Kernel functions

```@docs; canonical = false
compute_plant_available_water
```

```@docs; canonical = false
compute_plant_available_water!
```
