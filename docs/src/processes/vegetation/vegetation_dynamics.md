# Vegetation dynamics 

```@meta
CurrentModule = Terrarium
```

```@setup default
using Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Vegetation dynamics describes how the fractional area coverage $\nu \in [0, 1]$ of a plant functional type (PFT) within a grid cell changes over time. PFTs compete for the available space $(1 - \nu)$ on the landscape, with their ability to expand determined by productivity and limited by disturbance:

- **Productivity**: High net primary productivity (NPP) allows a PFT to expand into unoccupied space
- **Disturbance**: Fire, disease, and other mortality mechanisms reduce vegetation extent
- **Seed availability**: A minimum seed fraction ensures PFTs can reestablish after disturbance

Terrarium implements vegetation dynamics following the Lotka–Volterra approach of PALADYN (Willeit 2016), which balances growth and disturbance to determine PFT area fraction dynamics.

```@docs; canonical = false
AbstractVegetationDynamics
```

## PALADYN vegetation dynamics

```@docs; canonical = false
PALADYNVegetationDynamics
```

```@example default
variables(PALADYNVegetationDynamics(Float32))
```

### Lotka-Volterra vegetation competition

The vegetation area fraction ($\nu$) for a single PFT evolves according to:
```math
\begin{equation}
\frac{d\nu}{dt} = \left( \lambda_{\text{NPP}} \frac{C_{\text{veg}}}{\nu^*} \right) (1 - \nu) - \gamma_v \nu^*
\end{equation}
```

where:
- $\lambda_{\text{NPP}}$ is an NPP-dependent partitioning factor that controls how much productivity is available for expansion (see [Vegetation carbon dynamics](@ref))
- $C_{\text{veg}}$ is the vegetation carbon pool [kgC/m²]
- $\nu^* = \max(\nu, \nu_{\text{seed}})$ is the effective PFT fraction, bounded below by the seed fraction to prevent division by zero
- $\gamma_v$ is the disturbance rate [1/time]
- $(1 - \nu)$ is the fraction of the grid cell available for expansion

**Growth term**: $(\lambda_{\text{NPP}} C_{\text{veg}} / \nu^*) (1 - \nu)$ represents logistic growth, where productive PFTs with high carbon stocks expand rapidly into unoccupied space but slow as $\nu \to 1$.

**Disturbance term**: $\gamma_v \nu^*$ represents losses due to disturbance, proportional to PFT extent.

### Seed fraction and PFT establishment

The effective fraction $\nu^*$ is bounded below by the seed fraction $\nu_{\text{seed}}$ to ensure that a PFT at very low abundance can still persist and recover after disturbance, and to avoid division by zero in the growth term:
```math
\begin{equation}
\nu^* = \max(\nu, \nu_{\text{seed}})
\end{equation}
```

Typical values are $\nu_{\text{seed}} \approx 0.001$ (0.1%).

### Disturbance rates

The disturbance rate $\gamma_v$ aggregates losses from fire, pest outbreaks, disease, and senescence. A minimum baseline disturbance rate $\gamma_v^{\text{min}}$ is always enforced:
```math
\begin{equation}
\gamma_v \geq \gamma_v^{\text{min}}
\end{equation}
```

In the current implementation, $\gamma_v = \gamma_v^{\text{min}}$ (constant). Future implementations may couple disturbance to soil moisture, fire danger indices, or other vegetation state variables.

### Coupling to productivity

The growth term couples vegetation dynamics to the carbon cycle through $C_{\text{veg}}$ and $\lambda_{\text{NPP}}$, both computed by the carbon dynamics model (see [Vegetation carbon dynamics](@ref)). $\lambda_{\text{NPP}}$ acts as a gate: when a PFT's leaf area index is below its minimum viable threshold, $\lambda_{\text{NPP}} = 0$ and all NPP is directed toward building carbon rather than expanding area. Once the PFT is sufficiently productive, $\lambda_{\text{NPP}} > 0$ and area expansion becomes possible.

## Methods

```@docs; canonical = false
compute_γv
```

```@docs; canonical = false
compute_ν_star
```

```@docs; canonical = false
compute_ν_tendency
```

## Kernel functions

```@docs; canonical = false
compute_ν_tendencies!
```
