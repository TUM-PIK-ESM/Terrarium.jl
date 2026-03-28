# Skin temperature and ground heat flux

```@meta
CurrentModule = Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

The **skin temperature** (sometimes also called *surface temperature*) $T_s$ is the effective radiative temperature of the land surface. It represents the temperature at which the surface emits longwave radiation and is the temperature that directly controls evaporation and sensible heat flux in the surface energy balance.

Skin temperature differs from the air temperature at height $z$ by the effects of turbulent mixing and surface properties. It is not necessarily the same as subsurface soil temperature, $T_g$ (the temperature of the top soil layer), due to the thermal resistance between surface and soil.

## Implicit skin temperature

The implicit approach solves for skin temperature that simultaneously satisfies the surface energy balance:

```@docs; canonical = false
ImplicitSkinTemperature
```

```math
\begin{equation}
R_{\text{net}}(T_s) = H_s(T_s) + H_l(T_s) + G(T_s)
\end{equation}
```

This requires solving a nonlinear equation at each grid point and time step. The implicit approach is more physically consistent but computationally more intensive.

## Prescribed skin temperature

```@docs; canonical = false
PrescribedSkinTemperature
```

## Methods

```@docs; canonical = false
compute_skin_temperature
```

```@docs; canonical = false
update_skin_temperature!
```

```@docs; canonical = false
compute_ground_heat_flux!
```

```@docs; canonical = false
compute_ground_heat_flux
```
