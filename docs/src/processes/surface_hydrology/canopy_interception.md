# Canopy interception

```@meta
CurrentModule = Terrarium
```

```@setup interception
using Terrarium
using InteractiveUtils
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

When precipitation falls on vegetated ground, some fraction is intercepted and stored as liquid water on leaves and stems. This intercepted water either evaporates directly back to the atmosphere or eventually reaches the ground. Leaf absorption is also possible but typically constitutes a fairly small fraction in intercepted water. The fraction of falling precipitation intercepted by the canopy directly depends on vegetation density and structural properties (leaf area index and stem area index).

```@docs; canonical = false
AbstractCanopyInterception
```

```@example interception
subtypes(Terrarium.AbstractCanopyInterception)
```

## Canopy interception scheme from PALADYN

```@docs; canonical = false
PALADYNCanopyInterception
```

```@example interception
variables(PALADYNCanopyInterception(Float32))
```

### Interception fraction

Following PALADYN [willeitPALADYNV10Comprehensive2016](@cite), the fraction of precipitation intercepted is
```math
\begin{equation}
I_{\text{can}} = \alpha_{\text{int}} P (1 - e^{-k_{\text{ext}}(\text{LAI} + \text{SAI})})\,,
\end{equation}
```
where $\alpha_{\text{int}}$ is the interception factor (-), $P$ is incident precipitation (m/s), $k_{\text{ext}}$ is the radiation extinction coefficient (-), $\text{LAI}$ is the leaf area index (m²/m²), and $\text{SAI}$ is the stem area index (m²/m²).

The term $(1 - e^{-k_{\text{ext}}(\text{LAI} + \text{SAI})})$ represents the vegetation cover fraction following the Beer-Lambert law (see e.g. [vandijkModellingInterception2001](@cite)) and increases from 0 (bare ground) toward 1 (dense forest).

### Canopy water storage

The canopy water balance is expressed as
```math
\begin{equation}
\frac{\partial w_{\text{can}}}{\partial t} = I_{\text{can}} - E_{\text{can}} - R_{\text{can}}\,,
\end{equation}
```
where $w_{\text{can}}$ is liquid water stored on the canopy (m), $I_{\text{can}}$ is the interception rate (m/s), $E_{\text{can}}$ is evaporation of intercepted water (m/s), and $R_{\text{can}}$ is the removal rate (m/s).

The canopy storage capacity $w_{\text{can, max}}$ (m) is diagnosed from the sum of the leaf and stem area indices,
```math
\begin{equation}
w_{\text{can, max}} = w_0 (\text{LAI} + \text{SAI})\,,
\end{equation}
```
where $w_0$ (m) is the specific water capacity per unit leaf/stem area.

The removal rate of water from the canopy (e.g. due to gravity-induced drip) is computed as
```math
\begin{equation}
R_{\text{can}} = \frac{w_{\text{can}}}{\tau_w}\,,
\end{equation}
```
where $\tau_w$ is the canopy water removal timescale (typically 1 day = 86400 s) (s).

The saturation fraction of the canopy is
```math
\begin{equation}
f_{\text{can}} = \frac{w_{\text{can}}}{w_{\text{can, max}}}\,,
\end{equation}
```
which ranges from 0 (dry) to 1 (saturated) and controls the efficiency of canopy evaporation.

### Precipitation

After accounting for interception and canopy , the precipitation reaching the ground is
```math
\begin{equation}
P_{\text{ground}} = P - I_{\text{can}} + R_{\text{can}} - E_{\text{can}}\,.
\end{equation}
```

This represents the water available for soil infiltration and surface runoff.

## Process interface

```@docs; canonical = false
compute_auxiliary!(state, grid, canopy_interception::PALADYNCanopyInterception, atmos::AbstractAtmosphere)
```

```@docs; canonical = false
compute_tendencies!(state, grid, canopy_interception::PALADYNCanopyInterception, evapotranspiration::AbstractEvapotranspiration, args...)
```

## Methods

```@docs; canonical = false
canopy_water
```

```@docs; canonical = false
saturation_canopy_water
```

```@docs; canonical = false
rainfall_ground
```

## Kernel functions

```@docs; canonical = false
compute_canopy_interception
```

```@docs; canonical = false
compute_canopy_saturation_fraction
```

```@docs; canonical = false
compute_canopy_water_removal
```

```@docs; canonical = false
compute_w_can_tendency
```

```@docs; canonical = false
compute_precip_ground
```

## [References](@id "canopy_interception.refs")

```@bibliography
Pages = ["canopy_interception.md"]
Canonical = false
```
