# Snow parameterizations

```@meta
CurrentModule = Terrarium
```

```@setup snowparams
using Terrarium
using InteractiveUtils
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

[`SingleLayerSnow`](@ref) composes independently swappable sub-parameterizations for the areal coverage,
bulk density, thermal conductivity, and hydraulic properties of the pack. Together with the snow water
equivalent and bulk density these determine the geometric and thermal properties (snow depth, cover
fraction, thermal conductivity) diagnosed each auxiliary pass.

## Snow cover

The sub-grid snow-covered area fraction ``f_\text{snow}`` controls how strongly the snowpack modifies the
surface albedo, conduction, and latent flux, and how rainfall is partitioned between the pack and the
bare ground.

```@docs; canonical = false
FractionalSnowCover

compute_snow_cover_fraction
```

## Snow density

The bulk density ``\rho_s`` converts the snow water equivalent to a physical snow depth and sets the
thermal conductivity.

```@docs; canonical = false
ConstantSnowDensity

snow_density(snow::SingleLayerSnow)

compute_snow_depth
```

## Thermal conductivity

The bulk snow thermal conductivity is parameterized as a function of the bulk density. The default is the
power-law form of [yenReviewThermalProperties1981](@cite); logarithmic and piecewise-quadratic forms
following [sturmThermalConductivitySeasonal1997](@cite) are also provided.

```@docs; canonical = false
PowerLawSnowThermalConductivity

LogarithmicSnowThermalConductivity

QuadraticSnowThermalConductivity

compute_thermal_conductivity(snow::SingleLayerSnow, constants::MaterialConstants, ρ_s)
```

## Hydraulic properties

The hydraulic properties set the Darcy-type meltwater outflow from the pack (see [Snow mass balance](@ref)).

```@docs; canonical = false
ConstantSnowHydraulics
```

## Kernel functions

```@docs; canonical = false
compute_snow_properties!
```

## [References](@id "snowparams.refs")

```@bibliography
Pages = ["snow_parameterizations.md"]
Canonical = false
```
