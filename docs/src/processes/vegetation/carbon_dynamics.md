# Plant carbon dynamics and vegetation evolution

```@meta
CurrentModule = Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Theory

### Vegetation carbon pools

Plant biomass is partitioned among multiple carbon pools: leaves, stems, and roots. These pools have different turnover rates and play distinct roles in plant physiology and ecosystem function. The carbon dynamics system tracks the total vegetation carbon and derives the balanced leaf area index (LAI) assuming a fixed biomass allocation scheme.

### Leaf area index

The specific leaf area (SLA) and allometric relationships determine how vegetation carbon is distributed:
```math
\begin{equation}
\text{LAI}_b = \frac{C_{\text{leaf}}}{C_0} = \frac{1}{\text{SLA}} \left( C_{\text{veg}} - C_{\text{struct}} \right)
\end{equation}
```

where $C_{\text{veg}}$ is the total vegetation carbon, $C_{\text{struct}}$ is the structural carbon (stems and roots), and LAI$_b$ is the balanced leaf area index that would exist at equilibrium with the current carbon pool.

### Leaf turnover and litterfall

Leaf turn occurs at a constant annual rate $\gamma_L$ (year$^{-1}$), representing senescence and leaf drop. The litterfall rate integrating across the canopy is given by
```math
\begin{equation}
\Lambda = \gamma_L C_{\text{leaf}} = \gamma_L \text{SLA}^{-1} \left( C_{\text{veg}} - C_{\text{struct}} \right)
\end{equation}
```

Other carbon pools (roots and stems) turn over at their respective rates $\gamma_R$ and $\gamma_S$ (year$^{-1}$).

### Net primary productivity partitioning

Net primary production (NPP) is partitioned between two pathways: increasing vegetation carbon in existing vegetation and enabling the spread of new vegetation. The partitioning factor $\lambda_{\text{NPP}}$ depends on whether the current LAI is above the minimum LAI threshold:
```math
\begin{equation}
\lambda_{\text{NPP}} = \begin{cases}
0 & \text{if } \text{LAI}_b < \text{LAI}_{\min} \\
\frac{\text{LAI}_b - \text{LAI}_{\min}}{\text{LAI}_{\max} - \text{LAI}_{\min}} & \text{if } \text{LAI}_{\min} \leq \text{LAI}_b \leq \text{LAI}_{\max} \\
1 & \text{if } \text{LAI}_b > \text{LAI}_{\max}
\end{cases}
\end{equation}
```

This ensures that when vegetation is below its minimum viable LAI, all NPP goes to increasing carbon rather than spreading.

### Vegetation area fraction dynamics

The area fraction $\nu$ of a given plant functional type (PFT) evolves according to a logistic-growth style equation:
```math
\begin{equation}
\frac{\partial \nu}{\partial t} = \lambda_{\text{NPP}} \frac{C_{\text{veg}}}{\nu^*} (1 - \nu) - \gamma_v \nu^*
\end{equation}
```

where $\nu^* = \max(\nu, \nu_{\text{seed}})$ ensures that a minimum seed fraction is maintained, and $\gamma_v$ (year$^{-1}$) is the background disturbance/mortality rate that limits maximum vegetation coverage. The first term represents expansion driven by carbon accumulation, while the second represents loss from disturbance.

## Abstract types

```@docs; canonical = false
AbstractVegetationCarbonDynamics
```

```@docs; canonical = false
AbstractVegetationDynamics
```

## Concrete types

### Carbon Dynamics

```@docs; canonical = false
PALADYNCarbonDynamics
```

### Vegetation Dynamics

```@docs; canonical = false
PALADYNVegetationDynamics
```

## Methods

```@docs; canonical = false
compute_LAI_b
```

```@docs; canonical = false
compute_λ_NPP
```

```@docs; canonical = false
compute_litterfall_rate
```

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
compute_carbon_tendency
```

```@docs; canonical = false
compute_carbon_tendencies!
```
