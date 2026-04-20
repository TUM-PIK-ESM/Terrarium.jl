# Soil hydrology

```@meta
CurrentModule = Terrarium
```

```@setup soilwater
using Terrarium
using InteractiveUtils
```

## Overview

Soil hydrology processes characterize the dynamics of ground water in both saturated and unsaturated soil. It defines parameters and methods needed to compute water fluxes between layers and grid cells within the soil domain. Implementations should extend [`AbstractSoilHydrology`](@ref) and should generally consist of at least four components:
- A scheme for computing vertical water fluxes between soil layers
- A closure parameterization linking soil saturation and pressure head
- A parameterization for [`soil hydraulic properties`](@ref "Hydraulic properties")
- A forcing term representing user-defined, internal sources/sinks in the soil domain (not including evapotranspiration)

```@docs; canonical = false
AbstractSoilHydrology
```

```@example soilwater
subtypes(Terrarium.AbstractSoilHydrology)
```

Terrarium currently provides a single general implementation of `SoilHydrology` following the above interface:

```@docs; canonical = false
SoilHydrology
```

### [Process interface](@id soilhydrology.interface)

Dispatches for `NoFlow`:
```@docs; canonical = false
initialize!(state, grid, hydrology::SoilHydrology, soil::AbstractSoil, args...)

compute_auxiliary!(state, grid, hydrology::SoilHydrology, soil::AbstractSoil, args...)

compute_tendencies!(state, grid, hydrology::SoilHydrology, soil::AbstractSoil, args...)
```

Dispatches for `RichardsEq`:
```@docs; canonical = false
initialize!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    constants::PhysicalConstants,
    args...
) where {NF}

compute_auxiliary!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    args...
) where {NF}

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
variables(SoilHydrology(Float32, NoFlow()))
```

### Richardson-Richards equation for variably saturated flow

The vertical flow of water in porous media, such as soils, can be formulated as following the conservation law
```math
    \phi\frac{\partial \xi(\psi)}{\partial t} = - \boldsymbol{\nabla} \cdot \textbf{j}_{\text{w}} + F_{\text{w}}(z,t),
```
where (as defined in [Soil stratigraphy](@ref)) $\phi$ is the natural porosity (or saturated water content) of the soil volume and $\xi \in [0,1]$ the saturation of pore water/ice. $F_{\text{w}}(z,t)$ is an inhomogeneous source/sink (forcing) term (1/s) and $\textbf{j}_{\text{w}}$ the water flux vector (m/s). 

Vertical fluxes in the soil column can be represented by combining gravity-driven advection with Darcy's law
```math
\begin{equation}
\textbf{j}_{\text{w}}^{\mathrm{v}} =\textbf{j}_{\text{w}} \cdot \mathbf{n}_z = -\kappa_{\text{w}}\frac{\partial \left(\psi + z\right)}{\partial z},
\end{equation}
```
where $\psi$ is the matric potential (m), $\mathbf{n}_z$ the normal vector perpendicular to the surface and $\kappa_{\text{w}}$ the hydraulic conductivity (m/s). Given the positive upwards convention (see [Numerical core](@ref)), the $z-$axis is thus defined positive away from the surface. Substituting this equation into the aforementioned conservation law yields the widely known Richardson-Richards equation for variably saturated flow in porous media [richardsCapillaryConductionLiquids1931](@cite). When the divergence of the water flux vector ($\boldsymbol{\nabla} \cdot \textbf{j}_{\text{w}}$) is simplified to only consider vertical water fluxes, this gives:
```math
\begin{equation}
\phi\frac{\partial \xi(\psi)}{\partial t} = \frac{\partial}{\partial z}\left[\kappa_w \frac{\partial \left(\psi + z\right)}{\partial z}\right] + F_{\text{w}}(z,t)
\end{equation}
```
which is called the mixed-form equation, as the rate of change on the left-hand side is in $\xi$ but the vertical gradient on the right-hand side is in $\psi$ [bonanClimateChangeTEMSoilMoisture2019](@cite). 

```@docs; canonical = false
RichardsEq
```

```@example soilwater
variables(SoilHydrology(Float32, RichardsEq()))
```

## Hydraulic properties

Besides the static soil characteristics mentioned in [Soil stratigraphy](@ref) such as porosity ($\phi$), a key time-varying hydraulic property is the (non-saturated) hydraulic conductivity $\kappa_{\text{w}}$. Different ways of modelling $\kappa_{\text{w}}$ are defined as subyptes of [`AbstractUnsatK`](@ref). The simplest formulation, called [`UnsatKLinear`](@ref), models $\kappa_{\text{w}}$ as:
```math
\begin{equation}
\kappa_{\text{w}} = \frac{\theta_{\text{liq}}}{\theta + \theta_{\text{air}}} \kappa_{\text{w,sat}}\,
\end{equation}
```
with $\kappa_{\text{w,sat}}$ the saturated hydraulic conductivity (m/s).  A more complex formulation taking ice into account, called [`UnsatKVanGenuchten`](@ref),follows [vangenuchtenHydraulicConductivity1980](@cite), here presented in the form of [westermannCryoGridCommunityModel2023; Eq. (26)](@cite):
```math
\begin{equation}
\kappa_{\text{w}} = I_\text{ice}\left(\frac{\theta_{\text{liq}}}{\phi}\right)^{0.5} \left[1 - \left(1 - \left(\frac{\theta_{\text{liq}}}{\phi}\right)^{n/(n+1)}\right)^{(n-1)/n}\right]^2 \kappa_{\text{w,sat}}\,,
\end{equation}
```
which adds an ice impedance factor $I_\text{ice} = 10^{\Omega\theta_{\text{ice}}/\theta}$ to the classic van Genuchten formulation to account for the blocking of water-filled pores by ice. 

```@docs; canonical = false
UnsatKLinear
```

```@docs; canonical = false
UnsatKVanGenuchten
```

Values for $\kappa_{\text{w,sat}}$, wilting point ($\theta_{wp}$), and field capacity   ($\theta_{fc}$) can be set as constants (see [`ConstantSoilHydraulics`](@ref)) or can be calculated from soil texture using pedotransfer functions (PTFs). Currently, PTFs based on [noilhanISBA1996](@cite) for $\theta_{wp}$ and $\theta_{fc}$ are available in [`SoilHydraulicsSURFEX`](@ref).

```@docs; canonical = false
ConstantSoilHydraulics
```

```@docs; canonical = false
SoilHydraulicsSURFEX
```

## Closures

Critical to the modelling of water transport in the soil, is the relationship between soil water/ice saturation ($\xi$) and the matric potential ($\Psi_m$), which is called the soil-water retention curve (SWRC). Recall from [the documentation on the core interfaces](@ref core_interfaces_hydrology_closure) that the total hydraulic head $\Psi$  is defined as:
```math
\begin{equation}
\Psi = \Psi_m(\xi) + \Psi_z + \Psi_h
\end{equation}
```
with $\Psi_z$ the elevation head and $\Psi_h$ the hydrostatic head contributed by free water above the water table. Currently, two implementations of the SWRC are available in Terrarium via [FreezeCurves.jl](https://github.com/CryoGrid/FreezeCurves.jl): the Brooks-Corey ([`BrooksCorey`](@extref FreezeCurves.BrooksCorey)) and van Genuchten ([`VanGenuchten`](@extref FreezeCurves.VanGenuchten)) formulations. 

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
Pages = ["soil_hydrology.md"]
Canonical = false
```
