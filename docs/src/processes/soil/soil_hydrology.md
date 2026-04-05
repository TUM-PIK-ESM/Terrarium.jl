# Soil hydrology

```@meta
CurrentModule = Terrarium
```

```@setup soilwater
using Terrarium
```

## Overview

Soil hydrology processes characterize the dynamics of ground water in both saturated and unsaturated soil. It defines parameters and methods needed to compute water fluxes between layers and grid cells within the soil domain. Implementations should extend [`AbstractSoilHydrology`](@ref) and should generally consist of at least four components:
- A scheme for computing vertical water fluxes between soil layers
- A closure parameterization linking soil saturation and pressure head
- A parameterization for [`soil hydraulic properties`](@ref "Hydraulic properties")
- A forcing term representing user-defined, internal sources/sinks in the soil domain (not including evapotranspiration)

```@docs; canonical = false
SoilHydrology
```

### [Process interface](@id soilhydrology.interface)

```@docs; canonical = false
initialize!(state, grid, hydrology::SoilHydrology, soil::AbstractSoil, args...)
initialize!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    constants::PhysicalConstants,
    args...
) where {NF}
```

```@docs; canonical = false
compute_auxiliary!(state, grid, hydrology::SoilHydrology, soil::AbstractSoil, args...)
compute_auxiliary!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    args...
) where {NF}
```

```@docs; canonical = false
compute_tendencies!(state, grid, hydrology::SoilHydrology, soil::AbstractSoil, args...)
compute_tendencies!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    constants::PhysicalConstants,
    evapotranspiration::Optional{AbstractEvapotranspiration},
    runoff::Optional{AbstractSurfaceRunoff},
    args...
) where {NF}
```

## Vertical flow

### Static soil hydrology ("No Flow")

The simplest possible soil hydrology scheme is one in which the soil saturation state remains constant over time. This can be appropriate for simulations in regions where the soil is normally waterlogged or where changes in the hydrological state can be otherwise assumed to play a negligible role. 

```@docs; canonical = false
NoFlow
```

```@example soilwater
variables(SoilHydrology(Float32))
```

### Richardson-Richards equation for variably saturated flow

The vertical flow of water in porous media, such as soils, can be formulated as following the conservation law
```math
    \phi\frac{\partial\vartheta(\psi)}{\partial t} - \boldsymbol{\nabla} \cdot \textbf{j}_{\text{w}} - F_{\text{w}}(z,t) = 0,
```
where $\phi$ is the natural porosity (or saturated water content) of the soil volume and $F_{\text{w}}(z,t)$ is an inhomogeneous source/sink (forcing) term (m/s).

Vertical fluxes in the soil column be represented by combining gravity-driven advection with Darcy's law
```math
\begin{equation}
\textbf{j}_{\text{w}} \cdot \mathbf{n} = -\kappa_{\text{w}}\frac{\partial \left(\psi + z\right)}{\partial z},
\end{equation}
```
where $\psi$ is the matric potential (m). Substituting this equation into the aforementioned conservation law yields the widely known Richardson-Richards equation for variably saturated flow in porous media [richardsCapillaryConductionLiquids1931](@cite).

```@docs; canonical = false
RichardsEq
```

```@example soilwater
variables(SoilHydrology(Float32, RichardsEq()))
```

## Hydraulic properties

```@docs; canonical = false
ConstantSoilHydraulics
```

```@docs; canonical = false
SoilHydraulicsSURFEX
```

## Closures

```@docs; canonical = false
SoilSaturationPressureClosure
```

## Methods

```@docs; canonical = false
get_hydraulic_properties
```

```@docs; canonical = false
get_swrc
```

```@docs; canonical = false
get_closure
```

```@docs; canonical = false
compute_water_table!
```

```@docs; canonical = false
adjust_saturation_profile!
```

```@docs; canonical = false
compute_hydraulics!
```

## Kernel functions

```@docs; canonical = false
saturation_water_ice
```

```@docs; canonical = false
hydraulic_conductivity
```

```@docs; canonical = false
liquid_water_fraction
```

```@docs; canonical = false
water_table
```

```@docs; canonical = false
surface_excess_water
```

## [References](@id "soilhydrology.refs")

```@bibliography
Pages = ["processes/soil/soil_hydrology.md"]
Canonical = false
```
