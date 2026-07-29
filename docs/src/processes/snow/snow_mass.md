# Snow mass balance

```@meta
CurrentModule = Terrarium
```

```@setup snowmass
using Terrarium
using InteractiveUtils
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

The snow water equivalent ``W_\text{snow}`` (m) evolves according to a continuous mass balance,
```math
\begin{equation}
\frac{\partial W_\text{snow}}{\partial t} = P_s + R_\text{on snow} - M_r - E_\text{subl},
\end{equation}
```
where ``P_s`` is snowfall, ``R_\text{on snow} = f_\text{snow}\,R`` is the rainfall intercepted by the snow-covered fraction, ``M_r`` is the meltwater outflow draining from the snowpack base, and ``E_\text{subl}`` is the surface sublimation rate. In the coupled land model the bare-ground fraction ``(1 - f_\text{snow})`` of the rainfall reaches the soil directly and the meltwater outflow is added to it, forming the water input to the [Surface runoff](@ref) scheme.

## Meltwater outflow

Liquid water in excess of the capillary retention drains from the snowpack with a Darcy-type cubic
conductivity ([tarbotonSpatiallyDistributedEnergy1994](@cite), eqns 23–24, in excess-saturation form).

```@docs; canonical = false
compute_meltwater_outflow

snow_meltwater_flux
```

## Advected precipitation heat

Precipitation carries sensible and latent heat into the snowpack, referenced to liquid water at 0 °C (the
enthalpy reference of the closure, see [Snow energy balance](@ref)).

```@docs; canonical = false
compute_snow_precip_heat_flux
```

## Tendencies

```@docs; canonical = false
compute_snow_water_tendency

compute_snow_tendencies!
```

## [References](@id "snowmass.refs")

```@bibliography
Pages = ["snow_mass.md"]
Canonical = false
```
