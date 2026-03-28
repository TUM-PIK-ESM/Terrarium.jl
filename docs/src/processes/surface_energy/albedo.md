# Albedo and emissivity

```@meta
CurrentModule = Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Albedo is the fraction of incident shortwave radiation reflected back to space:

```math
\begin{equation}
\alpha = \frac{\text{SW}_{\uparrow}}{\text{SW}_{\downarrow}}
\end{equation}
```

where $S_{\uparrow}$ is the reflected (upwelling) shortwave and $\text{SW}_{\downarrow}$ is the incident (downwelling) shortwave. Albedo ranges from 0 (perfectly absorbing) to 1 (perfectly reflecting).

### Variability of albedo and emissivity

Albedo can depend on multiple factors:
- **Surface type**: Snow (0.7–0.9, highly reflective), vegetation (0.1–0.3), bare soil (0.2–0.4), water (0.05–0.15)
- **Soil moisture**: Darker, wetter soils have lower albedo
- **Vegetation density**: Denser vegetation lowers albedo
- **Snow age**: Fresh snow is bright; aging snow darkens with impurities

Emissivity is the efficiency with which the surface emits thermal (longwave) radiation:

```math
\begin{equation}
\epsilon = \frac{L_{\uparrow}}{\sigma T_0^4}
\end{equation}
```

where $L_{\uparrow}$ is the total upwelling longwave radiation (including both surface emission and reflection of incident longwave), and $\sigma T_0^4$ is the blackbody emission at skin temperature $T_0$.

Most natural surfaces behave fairly similarly to black bodies with emissivity ~0.9. However, snow and water covered surfaces tend to have higher emissivity (>0.95) while certain types of rock and manmade surfaces are known to have lower emissivities.

The total upwelling longwave includes both surface emission and reflection:
$$L_{\uparrow} = \epsilon \sigma T_0^4 + (1-\epsilon) L_{\downarrow}$$

## Implementations

Terrarium currently provides only two very simple formulations of albedo and emissivity. The simplest scheme is [`ConstantAlbedo`](@ref) which treats both albedo and emissivity as spatially and temporally constants parameters:

```@docs; canonical = false
ConstantAlbedo
```

Alternative, they can be prescribed as user-specified input fields via [`PrescribedAlbedo`](@ref):

```@docs; canonical = false
PrescribedAlbedo
```

In the case where they are prescribed, they are permitted to vary in time and/or space at the discretion of the user and/or the corresponding [`InputSource`](@ref). However, their values will not be derived from the model state and thus could lead to inconsistencies.

Dynamic (state-dependent) albedo and emissivity are not yet implemented but will be added soon!

## Methods

```@docs; canonical = false
albedo
```

```@docs; canonical = false
emissivity
```
