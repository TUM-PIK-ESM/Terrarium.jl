# Vegetation carbon dynamics

```@meta
CurrentModule = Terrarium
```

```@setup default
using Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Plant biomass is partitioned among multiple carbon pools: leaves, stems, and roots. These pools have different turnover rates and play distinct roles in plant physiology and ecosystem function. Carbon dynamics tracks the total vegetation carbon $C_{\text{veg}}$ [kgC/m²] as the single prognostic variable, with the balanced leaf area index LAI$_b$ computed from it as an auxiliary variable. The splitting of $C_{\text{veg}}$ into separate leaf, stem, and root pools is handled implicitly through fixed allometric relationships.

```@docs; canonical = false
AbstractVegetationCarbonDynamics
```

## PALADYN carbon dynamics

```@docs; canonical = false
PALADYNCarbonDynamics
```

```@example default
variables(PALADYNCarbonDynamics(Float32))
```

### Leaf area index

The balanced leaf area index LAI$_b$ represents the leaf area that would be in equilibrium with the current carbon pool. It is derived from $C_{\text{veg}}$ using the specific leaf area (SLA) and allometric relationships between leaf, stem, and root biomass (Eqs. 76–79, PALADYN, Willeit 2016):

```math
\begin{equation}
\text{LAI}_b = \frac{C_{\text{leaf}}}{C_0} = \frac{1}{\text{SLA}} \left( C_{\text{veg}} - C_{\text{struct}} \right)
\end{equation}
```

where $C_{\text{struct}}$ is the structural carbon (stems and roots) derived from $C_{\text{veg}}$ via allometric scaling, and $C_0$ is a normalisation constant. See [`compute_LAI_b`](@ref) for the implementation.

### Leaf turnover and litterfall

Leaf turnover occurs at a constant annual rate $\gamma_L$ (year$^{-1}$), representing senescence and leaf drop. Roots and stems turn over at their respective rates $\gamma_R$ and $\gamma_S$ (year$^{-1}$). The total local carbon loss rate $\Lambda_{\text{loc}}$, aggregating turnover across all pools, is (Eq. 75, PALADYN, Willeit 2016):

```math
\begin{equation}
\Lambda_{\text{loc}} = \left( \frac{\gamma_L}{\text{SLA}} + \frac{\gamma_R}{\text{SLA}} + \gamma_S \cdot a_{wl} \right) \text{LAI}_b
\end{equation}
```

where $a_{wl}$ is the allometric coefficient relating stem biomass to leaf area.

### Net primary productivity partitioning

NPP (computed by the autotrophic respiration model, see [Autotrophic respiration](@ref)) is partitioned between two pathways: increasing vegetation carbon in already-vegetated areas, and enabling PFT expansion into new area. The partitioning factor $\lambda_{\text{NPP}} \in [0, 1]$ depends on whether the current LAI$_b$ exceeds the minimum viable threshold (Eq. 74, PALADYN, Willeit 2016):

```math
\begin{equation}
\lambda_{\text{NPP}} = \begin{cases}
0 & \text{if } \text{LAI}_b < \text{LAI}_{\min} \\
\frac{\text{LAI}_b - \text{LAI}_{\min}}{\text{LAI}_{\max} - \text{LAI}_{\min}} & \text{if } \text{LAI}_{\min} \leq \text{LAI}_b \leq \text{LAI}_{\max} \\
1 & \text{if } \text{LAI}_b > \text{LAI}_{\max}
\end{cases}
\end{equation}
```

When $\text{LAI}_b < \text{LAI}_{\min}$, the PFT is below its minimum viable leaf area and all NPP is directed toward building carbon rather than expanding area. $\lambda_{\text{NPP}}$ is also consumed by the vegetation dynamics model to drive PFT area expansion (see [Vegetation dynamics](@ref)). The vegetation carbon tendency is then (Eq. 72, PALADYN, Willeit 2016):

```math
\begin{equation}
\frac{dC_{\text{veg}}}{dt} = (1 - \lambda_{\text{NPP}}) \cdot \text{NPP} - \Lambda_{\text{loc}}
\end{equation}
```

The first term is NPP retained for carbon accumulation on existing vegetated land; the second term represents losses through turnover of all pools.

## Methods

```@docs; canonical = false
compute_LAI_b
```

```@docs; canonical = false
compute_λ_NPP
```

```@docs; canonical = false
compute_Λ_loc
```

```@docs; canonical = false
compute_C_veg_tend
```

## Kernel functions

```@docs; canonical = false
compute_veg_carbon_tendency
```

```@docs; canonical = false
compute_veg_carbon_auxiliary!
```

```@docs; canonical = false
compute_veg_carbon_tendencies!
```

```@docs; canonical = false
compute_auxiliary_kernel!
```

```@docs; canonical = false
compute_tendencies_kernel!
```
