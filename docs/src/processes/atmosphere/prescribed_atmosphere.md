```@meta
CurrentModule = Terrarium
```

```@setup prescribed_atmosphere
using Terrarium
```

# Prescribed atmosphere

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

The atmosphere module in Terrarium provides the meteorological boundary conditions
for all surface energy, hydrological, and biogeochemical processes. A `PrescribedAtmosphere`
reads all forcing variables directly from input data — typically reanalysis products
(e.g. ERA5-Land) or observational records — rather than computing them interactively
with a coupled atmospheric model.

The primary state variables are:

| Variable | Symbol | Units |
|---|---|---|
| Near-surface air temperature | $T_a$ | °C |
| Atmospheric pressure | $p$ | Pa |
| Specific humidity | $q_a$ | kg/kg |
| Wind speed | $V_a$ | m/s |
| Rainfall | | m/s |
| Snowfall | | m/s |
| Downwelling shortwave radiation | $\text{SW}_{\downarrow}$ | W/m² |
| Downwelling longwave radiation | $\text{LW}_{\downarrow}$ | W/m² |
| Daytime length | | hr |
| CO₂ concentration | | ppm |

All of these are treated as spatially and temporally varying input fields that may be updated in
each timestep by their corresponding [`InputSource`](@ref)s.

### Vapor pressure deficit

The vapor pressure deficit (VPD) quantifies how far the atmosphere is from
saturation. It is used to drive evapotranspiration and is computed from the
air temperature, specific humidity, and atmospheric pressure:

```math
\begin{equation}
\Delta e = e_{\text{sat}}(T_s) - e_a(q_a, p)
\end{equation}
```

where $e_{\text{sat}}(T_s)$ is the saturation vapor pressure at surface temperature $T_s$,
and $e_a = q_a p / \varepsilon$ is the actual vapor pressure, with $\varepsilon \approx 0.622$
the ratio of molecular weights of water vapor to dry air. The corresponding specific humidity
deficit is

```math
\begin{equation}
\Delta q = \frac{\varepsilon \Delta e}{p}.
\end{equation}
```

### Aerodynamic resistance

The `PrescribedAtmosphere` delegates to an [`AbstractAerodynamics`](@ref) parameterization
to compute the aerodynamic resistance $r_a$ [s/m] between the land surface and the
atmosphere, with [`ConstantAerodynamics`](@ref) used by default. See [Aerodynamics](@ref) for
further discussion.

## Implementations

```@docs; canonical = false
PrescribedAtmosphere
```

`compute_auxiliary!` and `compute_tendencies!` are both no-ops for `PrescribedAtmosphere`:
all state variables are provided as inputs and updated by the input layer rather than
computed from process equations.

### Humidity parameterizations

```@docs; canonical = false
SpecificHumidity
```

### Precipitation parameterizations

```@docs; canonical = false
RainSnow
```

### Radiation parameterizations

```@docs; canonical = false
LongShortWaveRadiation
```

### Tracer gases

```@docs; canonical = false
TracerGas
```

```@docs; canonical = false
AmbientCO2
```

```@docs; canonical = false
TracerGases
```

## Abstract types

```@docs; canonical = false
AbstractAtmosphere
```

```@docs; canonical = false
AbstractHumidity
```

```@docs; canonical = false
AbstractPrecipitation
```

```@docs; canonical = false
AbstractIncomingRadiation
```

## Kernel functions

```@docs; canonical = false
air_temperature
```

```@docs; canonical = false
air_pressure
```

```@docs; canonical = false
windspeed
```

```@docs; canonical = false
specific_humidity
```

```@docs; canonical = false
rainfall
```

```@docs; canonical = false
snowfall
```

```@docs; canonical = false
shortwave_down
```

```@docs; canonical = false
longwave_down
```

```@docs; canonical = false
daytime_length
```

```@docs; canonical = false
aerodynamic_resistance
```

```@docs; canonical = false
compute_vpd
```

```@docs; canonical = false
compute_humidity_vpd
```
