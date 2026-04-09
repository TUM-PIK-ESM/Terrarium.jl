# Evapotranspiration

```@meta
CurrentModule = Terrarium
```

```@setup ET
using Terrarium
using InteractiveUtils
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Evapotranspiration (ET) is the combined process of water evaporation from soil and open water surfaces, evaporation of water intercepted by the canopy, and transpiration through leaf stomata. These processes remove water from the surface, driving latent heat flux and competing with sensible heat and ground heat fluxes in the surface energy balance.

All evapotranspiration pathways are primarily driven by the **vapor pressure deficit** or **specific humidity deficit** $\Delta q = q_{\text{sat}}(T_s) - q_a$, where $q_{\text{sat}}(T_s)$ is the saturation specific humidity at surface temperature $T_s$ (kg/kg) and $q_a$ is the atmospheric specific humidity at reference height (kg/kg).

Each pathway is also modulated by **aerodynamic resistance** $r_a$ (s/m) (between surface and atmosphere) and/or **stomatal resistance** $r_s$ (s/m) (in transpiration).

```@docs; canonical = false
AbstractEvapotranspiration
```

```@example ET
subtypes(Terrarium.AbstractEvapotranspiration)
```

## Bare ground evaporation

```@docs; canonical = false
BareGroundEvaporation
```

```@example ET
variables(BareGroundEvaporation(Float32))
```

## Canopy evapotranspiration

```@docs; canonical = false
PALADYNCanopyEvapotranspiration
```

```@example ET
variables(PALADYNCanopyEvapotranspiration(Float32))
```

### Evaporation from the canopy

Evaporation of water intercepted by the canopy depends on the saturation state of the canopy (fraction of leaves wet),
```math
\begin{equation}
E_{\text{can}} = f_{\text{can}} \frac{\Delta q}{r_a}\,,
\end{equation}
```
where $f_{\text{can}}$ is the canopy saturation fraction (0 = dry, 1 = saturated) (-) and $\Delta q$ is the vapor pressure deficit (kg/kg).

When $f_{\text{can}} = 0$ (completely dry canopy), $E_{\text{can}} = 0$. When $f_{\text{can}} = 1$ (wet canopy), evaporation proceeds at the potential rate.

### Ground evaporation

Evaporation from exposed soil or under-canopy surfaces is limited by soil water availability,
```math
\begin{equation}
E_{\text{ground}} = \beta \frac{\Delta q}{r_a + r_e}\,,
\end{equation}
```
where $\beta$ is the ground evaporation resistance factor (0 to 1) (-) and $r_e$ is the aerodynamic resistance between ground and canopy (s/m).

The resistance factor $\beta$ is computed from soil moisture in the upper layer: $\beta = 1$ when soil is wet (at field capacity) and $\beta \to 0$ as soil dries.

### Transpiration

Plant transpiration occurs through stomata and is controlled by stomatal conductance,
```math
\begin{equation}
T_{\text{canopy}} = \frac{\Delta q}{r_a + r_s}\,,
\end{equation}
```
where $r_s = 1 / g_w$ is the stomatal resistance (s/m) and $g_w$ is the stomatal conductance (m/s) (computed from photosynthesis; see ...).

High stomatal conductance (when photosynthetically active) leads to low stomatal resistance and high transpiration. This creates a strong coupling between carbon uptake (photosynthesis) and water loss (transpiration).

### Total evapotranspiration

The PALADYN approach combines all three pathways,
```math
\begin{equation}
\text{ET} = E_{\text{can}} + E_{\text{ground}} + T_{\text{canopy}}\,,
\end{equation}
```
to obtain a total surface humidity flux that can be converted into the [latent heat flux][@ref "Latent heat flux"] for use in the [surface energy balance](@ref "Surface energy balance").

## Process interface

```@docs; canonical = false
compute_auxiliary!(state, grid, ::BareGroundEvaporation, ::NoCanopyInterception, ::PhysicalConstants, ::AbstractAtmosphere, ::Optional{AbstractSoil}, args...)
```

```@docs; canonical = false
compute_auxiliary!(state, grid, ::PALADYNCanopyEvapotranspiration, ::AbstractCanopyInterception, ::PhysicalConstants, ::AbstractAtmosphere, ::AbstractSoil, ::AbstractVegetation, args...)
```

## Ground resistance parameterizations

```@docs; canonical = false
ConstantEvaporationResistanceFactor
```

```@docs; canonical = false
SoilMoistureResistanceFactor
```

## Coupling to soil hydrology

Subtypes of `AbstractEvapotranpsiration` automatically inherit an implementation of the [`forcing`](@ref) interface for [`SoilHydrology`](@ref),
which computes the contribution of ET to the soil moisture tendency in each soil layer. The default implementation draws the rescaled surface humidity
flux only from the uppermost soil layer.

```@docs; canonical = false
forcing(i, j, k, grid, clock, fields, evapotranspiration::AbstractEvapotranspiration, ::AbstractSoilHydrology)
```

## Kernel functions

```@docs; canonical = false
surface_humidity_flux
```

```@docs; canonical = false
compute_transpiration
```

```@docs; canonical = false
compute_evaporation_ground
```

```@docs; canonical = false
compute_evaporation_canopy
```

```@docs; canonical = false
ground_evaporation_resistance_factor
```
