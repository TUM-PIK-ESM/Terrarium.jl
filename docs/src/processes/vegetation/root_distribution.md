# Root distribution

```@meta
CurrentModule = Terrarium
```

```@setup vegroots
using Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

The vertical distribution of plant roots determines the soil profile from which a plant extracts water and nutrients. Root distributions are typically concentrated near the surface (where nutrient availability is highest) but extend to depth to access water during drier periods.

The root distribution can be described as a density function over the vertical (elevation) coordinate $z$ (m),
```math
\begin{equation}
r(z) = \frac{\partial R}{\partial z}
\end{equation}
```
where $R(z)$ is the cumulative fraction of roots from the surface to $z$. The total roots integrate to unity: $\int_0^{z_{\max}} r(z) \, dz = 1$. Since the soil layers 

```@docs; canonical = false
AbstractRootDistribution
```

### Static exponential root distribution

A common choice for $r(z)$ is the average of two exponential distributions with different depth scales following [zengGlobalVegetationRoot2001](@cite),
```math
\begin{equation}
r(z) = \frac{1}{2} \left[ a \exp(a z) + b \exp(b z) \right]
\end{equation}
```
where $a$ and $b$ are empirical rate parameters controlling the depth e-folding scales (m⁻¹). The first parameter $a$ can be viewed as controlling shallow roots, while $b$ controls deeper root penetration. Larger values produce shallower distributions. Note that, in this formulation, we assume `z` to be *decreasing* with depth (negative below the surface).

Root distributions vary among plant functional types. Shallow-rooted plants (grasses, shrubs) may extract water almost entirely from the top soil layers, while deep-rooted trees access water at greater depths. These differences are captured through PFT-specific values of the rate parameters $a$ and $b$.

```@docs; canonical = false
StaticExponentialRootDistribution
```

```@example vegroots
variables(StaticExponentialRootDistribution(Float32))
```

## Methods

```@docs; canonical = false
root_density
```

```@docs; canonical = false
root_fraction
```

## [References](@id "rootdist.refs")

```@bibliography
Pages = ["root_distribution.md"]
Canonical = false
```
