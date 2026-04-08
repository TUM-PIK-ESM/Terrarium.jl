```@meta
CurrentModule = Terrarium
```

```@setup aerodynamics
using Terrarium
using InteractiveUtils
```

# Aerodynamics

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Aerodynamic processes govern the turbulent exchange of heat, water vapor, and momentum between the land surface and the overlying atmosphere. Surface fluxes are often approximated by assuming that the turbulent eddies act like a diffusive process with an effective resistance $r_a$ (s/m),
```math
\begin{equation}
\text{Flux} = \frac{\Delta X}{r_a}
\end{equation}
```
where $\Delta X$ is the gradient of the transported quantity (temperature, specific humidity, etc.) between the surface and the reference height. $r_a$ is typically referred to as *aerodynamic resistance* and can be equivalently formulated as a conductance $g_a = r_a^{-1}$ (s/m). This resistance term constitutes a key constraint on energy and water exchange between the atmosphere and the land surface. It can be related to the near-surface wind speed $V_a$ via the bulk drag coefficient $C_h$,
```math
\begin{equation}
r_a = \frac{1}{C_h V_a}.
\end{equation}
```
The drag coefficient $C_h$ absorbs the effects of surface roughness, the height of the reference level, and in more sophisticated schemes, the atmospheric stability of the boundary layer. In neutral conditions over smooth surfaces, typical values are $C_h \approx 10^{-3}$ to
$5 \times 10^{-3}$.

In Terrarium, wind speed is typically clipped to a small positive minimum $V_{\min}$ (see [`PrescribedAtmosphere`](@ref)) to prevent division by zero in calm conditions.

## Constant aerodynamics

In the simplest case, we can treat $C_h$ as a parameter that is invariant in space and time. This typically does not result in very realistic predictions of surface energy and water fluxes but can be a suitable approximation for idealized simulations and testing.

```@docs; canonical = false
ConstantAerodynamics
```

## Kernel functions

```@docs; canonical = false
drag_coefficient
```
