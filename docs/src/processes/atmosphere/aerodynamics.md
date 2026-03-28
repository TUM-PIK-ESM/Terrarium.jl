```@meta
CurrentModule = Terrarium
```

```@setup aerodynamics
using Terrarium
```

# Aerodynamics

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Aerodynamic processes govern the turbulent exchange of heat, water vapor, and momentum
between the land surface and the overlying atmosphere. In the bulk aerodynamic
formulation used by Terrarium, surface fluxes are approximated by assuming that the
turbulent eddies act like a diffusive process with an effective resistance $r_a$ [s/m]:

```math
\begin{equation}
\text{Flux} = \frac{\Delta X}{r_a}
\end{equation}
```

where $\Delta X$ is the gradient of the transported quantity (temperature, specific
humidity, etc.) between the surface and the reference height. This resistance
is related to the bulk drag coefficient $C_h$ and the near-surface wind speed $V_a$ by:

```math
\begin{equation}
r_a = \frac{1}{C_h V_a}
\end{equation}
```

The drag coefficient $C_h$ absorbs the effects of surface roughness, the height of the
reference level, and — in more sophisticated schemes — the atmospheric stability. In
neutral conditions over smooth surfaces, typical values are $C_h \approx 10^{-3}$ to
$5 \times 10^{-3}$.

### Aerodynamic conductance

The aerodynamic conductance $g_a = 1 / r_a = C_h V_a$ [m/s] is the inverse of the
resistance and often appears in flux expressions alongside stomatal and soil surface
conductances:

```math
\begin{equation}
g_a = C_h V_a
\end{equation}
```

Wind speed is clipped to a small positive minimum $V_{\min}$ (default 0.01 m/s) in
`PrescribedAtmosphere` to prevent division by zero in near-calm conditions.

## Implementations

```@docs; canonical = false
ConstantAerodynamics
```

`compute_auxiliary!` and `compute_tendencies!` are not defined for
`AbstractAerodynamics` implementations — aerodynamic quantities are computed on demand
via [`drag_coefficient`](@ref) and [`aerodynamic_resistance`](@ref).

## Kernel functions

```@docs; canonical = false
drag_coefficient
```
