# Radiative fluxes

```@meta
CurrentModule = Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

The radiative budget is characterized by the incidental downwelling (incoming) and upwelling (outgoing) radiative fluxes, each of which is defined in Terrarium as strictly nonnegative in their respective directions. The net radiative flux can then be computed as the sum of the directional fluxes for both the shortwave (solar) and longwave (thermal) components:

```math
\begin{equation}
R_{\text{net}}(T_s) = \text{SW}_{\uparrow} - \text{SW}_{\downarrow} + \text{LW}_{\uparrow}(T_s) - \text{LW}_{\downarrow}
\end{equation}
```

where $\text{SW}_{\uparrow} = \alpha \text{SW}_{\downarrow}$ is upwelling (reflected) shortwave radiation (W/m²), $\text{SW}_{\downarrow}$ is downwelling (incident) shortwave radiation (W/m²), $\text{LW}_{\uparrow} = \epsilon \sigma T_0^4 + (1-\epsilon) L_{\downarrow}$ is upwelling longwave radiation from the surface and reflected downwelling longwave (W/m²), $\text{LW}_{\downarrow}$ is downwelling (incident) longwave radiation (W/m²), $\epsilon$ is the surface emissivity (see [Albedo and Emissivity](albedo.md)) (-), $\alpha$ is the surface albedo (see [Albedo and Emissivity](albedo.md)) (-), $\sigma$ is the Stefan-Boltzmann constant (W/m²/K⁴), and $T_s$ is the skin temperature (K).

### Shortwave radiation

Shortwave radiation is determined by the solar angle, atmospheric transmittance, and cloud cover, all of which are assumed to be implicitly accounted for by the forcing data or a coupled atmosphere model. The absorbed shortwave depends on surface albedo $(1 - \alpha)$, which may vary with surface type (soil, vegetation, snow) and solar zenith angle.

### Longwave radiation

The surface emits longwave radiation according to the Stefan-Boltzmann law, with the emission temperature determined by the skin temperature $T_0$. The amount of downwelling longwave radiation is controlled by atmospheric water vapor, cloud cover, air temperature, and greenhouse gases.

## Diagnosed radiative fluxes

The standard implementation of the radiative balance is provided through [`DiagnosedRadiativeFluxes`](@ref):

```@docs; canonical = false
DiagnosedRadiativeFluxes
```

which computes the radiative budget according to the above equations.

## Prescribed radiative fluxes

Alternatively, the outgoing shortwave and longwave radiation can also be prescribed as input fields via

```@docs; canonical = false
PrescribedRadiativeFluxes
```

This may be useful for testing or in cases where you want to couple Terrarium with another model that will take care of computing the incoming and outgoing radiative fluxes.

## Process interface

```@docs; canonical = false
compute_auxiliary!(state, grid, rad::PrescribedRadiativeFluxes, seb::AbstractSurfaceEnergyBalance, atmos::AbstractAtmosphere, args...)
```

```@docs; canonical = false
compute_auxiliary!(state, grid, rad::DiagnosedRadiativeFluxes, seb::AbstractSurfaceEnergyBalance, consts::PhysicalConstants, atmos::AbstractAtmosphere, args...)
```

## Methods

```@docs; canonical = false
compute_shortwave_up
```

```@docs; canonical = false
compute_longwave_up
```

```@docs; canonical = false
compute_surface_net_radiation
```

## Kernel functions

```@docs; canonical = false
compute_surface_upwelling_radiation
```

```@docs; canonical = false
compute_radiative_fluxes!
```
