# Phenology

```@meta
CurrentModule = Terrarium
```

```@setup default
using Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Phenology describes the seasonal emergence and senescence of leaves. Models typically distinguish between three phenology types: evergreen, in which plants maintain constant leaf foliage throughout the year, summergreen (seasonal deciduous), in which leaves are present during the warm season and drop in the cold season, and raingreen (stress deciduous), in which leaves are present during the rainy season and drop in the dry season.


```@docs; canonical = false
AbstractPhenology
```

## PALADYN phenology model

```@docs; canonical = false
PALADYNPhenology
```

```@example default
variables(PALADYNPhenology(Float32))
```

This implementation follows the phenology scheme of PALADYN (Willeit, 2016), in which raingreen phenology is not represented. The phenology factor $\phi$ (0 to 1) represents the current fraction of the maximum leaf coverage, and $f_{\text{deciduous}}$ (0 to 1) is a climate-dependent smooth transition parameter between evergreen and deciduous behavior.


### Leaf area index computation

The leaf area index (LAI) is computed from the balanced LAI $\text{LAI}_b$ as follows

```math
\begin{equation}
\text{LAI} = (f_{\text{deciduous}} \cdot \phi + (1 - f_{\text{deciduous}})) \cdot \text{LAI}_b
\end{equation}
```

!!! warning
    Phenology is not fully implemented yet: currently $\phi = 1$ and $f_{\text{deciduous}} = 0$ which assumes an evergreen phenology.

## Methods

```@docs; canonical = false
compute_f_deciduous
```

```@docs; canonical = false
compute_phenology_factor
```

```@docs; canonical = false
compute_LAI
```

## Kernel functions

```@docs; canonical = false
compute_phenology
```

```@docs; canonical = false
compute_phenology!
```
