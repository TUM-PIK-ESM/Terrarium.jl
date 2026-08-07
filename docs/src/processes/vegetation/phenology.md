# Phenology

```@meta
CurrentModule = Terrarium
```

```@setup vegphen
using Terrarium
using InteractiveUtils
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Phenology describes the seasonal emergence and senescence of leaves. Models typically distinguish between three phenology types: evergreen, in which plants maintain constant leaf foliage throughout the year, summergreen (seasonal deciduous), in which leaves are present during the warm season and drop in the cold season, and raingreen (stress deciduous), in which leaves are present during the rainy season and drop in the dry season.


```@docs; canonical = false
AbstractPhenology
```

```@example vegphen
subtypes(Terrarium.AbstractPhenology)
```

## PALADYN phenology model

```@docs; canonical = false
PALADYNPhenology
```

```@example vegphen
variables(PALADYNPhenology(Float32))
```

This implementation follows the phenology scheme of PALADYN [willeitPALADYNV10Comprehensive2016](@cite), in which raingreen phenology is not represented. The phenology factor $\phi$ represents the current fraction of the maximum leaf coverage (0 to 1), and $f_{\text{deciduous}}$ is a climate-dependent smooth transition parameter (0 to 1) between evergreen and deciduous behavior.


### Leaf area index computation

The leaf area index (LAI) is computed from the balanced LAI $\text{LAI}_b$ as follows

```math
\begin{equation}
\text{LAI} = (f_{\text{deciduous}} \cdot \phi + (1 - f_{\text{deciduous}})) \cdot \text{LAI}_b
\end{equation}
```

!!! warning
    Phenology is not fully implemented yet: currently $\phi = 1$ and $f_{\text{deciduous}} = 0$ which assumes an evergreen phenology.

## Prescribed phenology

```@docs; canonical = false
PrescribedPhenology
```

```@example vegphen
variables(PrescribedPhenology(Float32))
```

The [`PrescribedPhenology`](@ref) scheme treats the leaf area index as an externally imposed (and possibly time-varying) input rather than deriving it from a prognostic carbon pool. It is intended for simulations driven by observed or reanalysis LAI (e.g. an ERA5-Land LAI climatology).

When plant traits are supplied by the enclosing vegetation component, the scheme additionally derives the phenology factor $\phi$ by inverting the deciduous relation $\text{LAI} = \phi \cdot \text{LAI}_b$, using the maximum (reference, annual-maximum) leaf area index $\text{LAI}_{\max}$ as the balanced LAI:

```math
\begin{equation}
\phi = \mathrm{clamp}\!\left(\frac{\text{LAI}}{\text{LAI}_{\max}},\ 0,\ 1\right)
\end{equation}
```

so that `PrescribedPhenology` becomes a drop-in producer of `phenology_factor` for downstream processes. $\text{LAI}_{\max}$ is a plant functional-type parameter carried by the vegetation traits.

!!! note
    The seasonal amplitude of $\phi$ is set entirely by how well $\text{LAI}_{\max}$ matches the true annual maximum of the prescribed LAI series. Because $\phi$ here is a structural ratio rather than a phenological leaf-out state, it is most appropriate for deciduous vegetation; for an evergreen plant functional type $\phi$ should remain near 1 regardless of LAI.

## Process interface

```@docs; canonical = false
compute_auxiliary!(state, grid, phenol::PALADYNPhenology, vegcarbon::PALADYNCarbonDynamics, atmos::AbstractAtmosphere)
```

## Methods

```@docs; canonical = false
compute_phenology_factor
```

```@docs; canonical = false
compute_leaf_area_index
```

## Kernel functions

```@docs; canonical = false
compute_phenology
```

```@docs; canonical = false
compute_phenology!
```

## [References](@id "phenology.refs")

```@bibliography
Pages = ["phenology.md"]
Canonical = false
```