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

The root distribution can be described as a probability density function over depth $z$:
```math
\begin{equation}
r(z) = \frac{\partial R}{\partial z}
\end{equation}
```

where $R(z)$ is the cumulative fraction of roots from the surface to depth $z$. The total roots integrate to unity: $\int_0^{z_{\max}} r(z) \, dz = 1$.

```@docs; canonical = false
AbstractRootDistribution
```

### Coupling to water availability

The root distribution is used to integrate soil water availability across the rooting profile:
```math
\begin{equation}
\beta(t) = \int_0^{z_{\max}} W(z, t) \cdot r(z) \, dz
\end{equation}
```

This provides the soil moisture limitation factor $\beta$ that constrains photosynthesis and transpiration.

### Static exponential root distribution

```@docs; canonical = false
StaticExponentialRootDistribution
```

```@example vegroots
variables(StaticExponentialRootDistribution(Float32))
```

A common parameterization based on observations is an average of two exponential distributions with different depth scales:
```math
\begin{equation}
r(z) = \frac{1}{2} \left[ a \exp(a z) + b \exp(b z) \right]
\end{equation}
```

where $a$ and $b$ are empirical rate parameters controlling the depth e-folding scales (m⁻¹). The first parameter $a$ typically controls shallow roots, while $b$ controls deeper root penetration. Larger values produce shallower distributions.

Root distributions vary among plant functional types. Shallow-rooted plants (grasses, shrubs) may extract water almost entirely from the top soil layers, while deep-rooted trees access water at greater depths. These differences are captured through PFT-specific values of the rate parameters $a$ and $b$.

## Methods

```@docs; canonical = false
root_density
```

```@docs; canonical = false
root_fraction
```
