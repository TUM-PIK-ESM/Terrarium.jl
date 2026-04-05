# Stomatal conductance

```@meta
CurrentModule = Terrarium
```

```@setup vegstomcond
using Terrarium
```


!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Stomata regulate gas exchange between the leaf and atmosphere, directly controlling the trade-off between carbon uptake during photosynthesis and water loss through transpiration. Similar to photosynthesis, stomatal conductance, $g_w$, depends on multiple environmental factors, including light availability, CO₂ concentration, and temperature. 

```@docs; canonical = false
AbstractStomatalConductance
```

### The Medlyn stomatal conductance model

```@docs; canonical = false
MedlynStomatalConductance
```

```@example vegstomcond
variables(MedlynStomatalConductance(Float32))
```

This implementation uses the optimal stomatal conductance model of Medlyn (2011) [medlynReconcilingOptimalEmpirical2011](@cite), adapted from PALADYN [willeitPALADYNV10Comprehensive2016](@cite), which derives stomatal conductance from water-use efficiency optimization as follows

```math
\begin{equation}
g_w = g_0 + 1.6 \frac{A_n}{c_a}  \left(1 + \frac{g_1}{\sqrt{\text{VPD}}}\right) 
\end{equation}
```

where $g_0$ is the minimum stomatal conductance,  $g_1$ is a PFT-specific slope parameter, VPD is the vapor pressure deficit, $A_n$ is the net photosynthesis and $c_a$ is the atmospheric CO₂ concentration.

$g_w$ and $A_n$ are also
related by the diffusion equation

```math
\begin{equation}
g_w = g_0 + 1.6 \frac{A_n}{c_a - c_i} 
\end{equation}
```

where $c_i$ is the intercellular CO2 concentration.

The ratio of intercellular to atmosphere CO₂ concentration $\lambda_c$ can then be derived as

```math
\begin{equation}
\lambda_c = 1 - \frac{1}{1 + \frac{g_1}{\sqrt{\text{VPD}}}}
\end{equation}
```

## Methods

```@docs; canonical = false
compute_gw_can
```

```@docs; canonical = false
compute_λc
```

## Kernel functions

```@docs; canonical = false
compute_stomatal_conductance!
```

```@docs; canonical = false
compute_stomatal_conductance
```

## [References](@id "stomatal_conductance.refs")

```@bibliography
Pages = ["processes/vegetation/stomatal_conductance.md"]
Canonical = false
```
