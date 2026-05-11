# Soil energy balance

```@meta
CurrentModule = Terrarium
```

```@setup soilenergy
using Terrarium
using InteractiveUtils
```

## Overview

Heat transfer in the subsurface can be represented according to the heat equation, with the upper boundary set to surface temperature and the lower boundary set to a constant positive heat flux representing heat produced by the inner earth [jaegerApplicationTheoryHeat1965,lachenbruchChangingClimateGeothermal1986](@cite). If both the upper and lower boundaries are assumed to be constant over time, the steady-state temperature profile takes the form of a continuous piecewise linear function increasing over depth with the slope determined by the thermal properties of the ground material. The instantaneous temperature field can then be generally represented as
```math
\begin{equation}
T(z,t) = T_0 + \frac{Q_{\text{geo}}}{\kappa_{\text{h}}(z)}z + \Delta T(z,t)
\end{equation}
```
where $T(z,t)$ is the temperature field (K) over depth $z$ (m) and time $t$ (s), $T_0$ is the mean annual GST (K), $Q_{\text{geo}}$ is the geothermal heat flux (W/m²), and $\kappa_{\text{h}}(z)$ is the thermal conductivity which may vary with depth depending on the material (W/m K). The last term $\Delta T(z,t)$ represents transient disturbances to the steady state temperature profile due to both seasonal and long-term fluctuations in the upper and lower boundary conditions of the vertical domain. Simulating the impacts of these transient changes is one of the primary objectives of most numerical permafrost and land surface models.

Diffusive heat flow in a solid medium is governed by Fourier's law,
```math
\begin{equation}
     \mathbf{j}_\text{h}^{\mathrm{v}} = \mathbf{j}_\text{h} \cdot \mathbf{n}_z = -\kappa_{\text{h}}\frac{\partial T}{\partial z}\,,
\end{equation}
```
where $\mathbf{j}_\text{h}$ is the diffusive heat flux vector (W/m²) and $\mathbf{n}_z$ is the upward facing normal vector along the vertical $z$ axis.


### Phase change of pore water/ice

Since ground materials are often porous, i.e., there exists void space between the solid particles, it is necessary to consider the potential presence of water and/or ice in this void space, which is hereafter referred to as pore space, or simply, soil pores. The thermal effects of water and ice can be accounted for by considering not only the temperature of the material but rather the total internal energy of the elementary volume. Combining the diffusive flux with a potential advective heat flux $\mathbf{j}_{\text{h}}^{\text{w}}$ due to water flow yields the energy conservation law,
```math
\begin{equation}
    \label{eq:energyconservation}
\frac{\partial U(T,\theta)}{\partial t} = - \boldsymbol{\nabla} \cdot \left(\mathbf{j}_\text{h} + \mathbf{j}_h^{\text{w}}\right) + F_h(z,t),
\end{equation}
```
where $U(T,\theta)$ (J/m³) is the volumetric internal energy as a function of temperature and total water/ice content $\theta$ (m³/m³), and $F_h(z,t)$ is an inhomogeneous heat source/sink (forcing) term.

The advective heat flux $\mathbf{j}_{\text{h}}^{\text{w}}$ can be represented as,
```math
\begin{equation}
\mathbf{j}_{\text{h}}^{\text{w}} = \left( c_{\text{w}} T + L_{\text{sl}} \right) \mathbf{j}_{\text{w}} \rho_{\text{w}}
\end{equation}
```
where $L_{\text{sl}}$ (J/kg) and $c_{\text{w}}$ (J/(kg K)) represent the specific latent heat of fusion (solid → liquid) and heat capacity of liquid water respectively. Note that $\mathbf{j}_{\mathrm{w}}$ is the water flux vector as defined in [Soil hydrology](@ref). This flux term accounts for the energy transferred by the movement of water within the soil matrix. In model configurations that neglect subsurface water flow, this flux term is implicitly assumed to be zero.

!!! todo "Advective heat flux"
    The current implementation does not yet consider the advective heat flux, but this will be added soon!

```@docs; canonical = false
SoilEnergyBalance
```

As for [Soil hydrology](@ref), only the vertical fluxes are currently considered. This simplifies the notation of ``\eqref{eq:energyconservation}`` to
```math
\begin{equation}
\frac{\partial U(T,\theta)}{\partial t} = \frac{\partial}{\partial z}\left[\kappa_{\text{h}} \frac{\partial T}{\partial z}\right] + F_h(z,t)\,.
\end{equation}
```

## [Process interface](@id soilenergy.dispatches)

```@docs; canonical = false
initialize!(state, grid, energy::SoilEnergyBalance, soil::AbstractSoil, constants::PhysicalConstants, args...)

compute_auxiliary!(state, grid, energy::SoilEnergyBalance, soil::AbstractSoil, args...)

compute_tendencies!(state, grid, energy::SoilEnergyBalance, soil::AbstractSoil, args...)
```

## Closures

The constitutive relationship between energy and temperature plays a critical role in simulating the thermal state of porous media. This relation can be defined in integral form as
```math
\begin{equation}
    U(T,\theta) = \int_{T_{\text{ref}}}^T \tilde{C}(x,\theta) \, \mathrm{d}x 
\end{equation}
```
where $\tilde{C}$ is referred to as the *effective* or *apparent* heat capacity and $T_{\text{ref}}$ is a reference temperature. Alternatively, the internal energy is defined following [dallamicoEnergyConservingFreezingSoil2011; Eq. (4)](@cite) as:
```math
\begin{equation}
U(T,\theta) = C(\theta_{\text{w}},\theta) (T - T_{\text{ref}}) + \rho_{\text{w}} L_{\text{sl}} \theta_{\text{w}}(\theta, T),
\end{equation}
```
where $\theta_{\text{w}}(T,\theta)$ is the volumetric unfrozen water content as a function of temperature $T$ and total water/ice content $\theta$; $C(\theta_{\text{w}},\theta)$  (J/(K m³)) is the bulk volumetric material heat capacity of the volume as a function of the unfrozen and total water contents;  $\rho_{\text{w}}$ corresponds to the density (kg/m³) of water. The apparent heat capacity is then defined as the derivative of the energy-temperature relation,
```math
\begin{equation}
\tilde{C}(T,\theta) := \frac{\partial U}{\partial T} =
\overbrace{C(\theta_{\text{w}},\theta) + T \frac{\partial C}{\partial \theta_{\text{w}}}\frac{\partial \theta_{\text{w}}}{\partial T}}^{\text{Sensible}} \,+\,
\overbrace{\rho_{\text{w}} L_{\text{sl}} \frac{\partial\theta_{\text{w}}}{\partial T}}^{\text{Latent}}\,,
\end{equation}
```
where the chain-rule is applied on $C(\theta_{\text{w}},\theta) (T - T_{\text{ref}})$ as both $C$ and $\theta_{\text{w}}$ are functions of $T$. The grouping of terms on the right-hand side show the partitioning of energy change into **sensible** and **latent** heat. The sensible component represents the energy necessary to heat a volume of the material to a particular temperature, whereas the latent component corresponds to the energy required for the phase change of water in the volume from solid (frozen) to liquid (thawed).

In the simplest case where we neglect the effect of capillary action in the soil, the energy-temperature relation can be derived according to that of "free" water (i.e. unbound by the soil matrix),
```math
\begin{equation}
    \theta_{\text{w}}(U) =
        \begin{cases}
            0                   & U < -\rho_{\text{w}}L_{\text{sl}}\theta \\
            \frac{U}{L} & -\rho_{\text{w}}L_{\text{sl}}\theta \leq U < 0 \\
            \theta              & U \geq 0\,,
        \end{cases}
\end{equation}
```
with temperature then determined by
```math
\begin{equation}
    T = U^{-1}(U(T,\theta)) =
    \begin{cases}
    \frac{U(T,\theta) - \rho_{\text{w}}L_{\text{sl}}\theta}{C} & U(T,\theta) < -\rho_{\text{w}}L_{\text{sl}}\theta \\
    0 & 0 \leq U(T,\theta) \leq \rho_{\text{w}}L_{\text{sl}}\theta \\
    \frac{U(T,\theta)}{C} &   U(T,\theta) \geq 0\,,
    \end{cases}.
\end{equation}
```

For more information on more complex formulations of these so-called soil freezing characteristic curves, see [FreezeCurves.jl](https://github.com/CryoGrid/FreezeCurves.jl/tree/main?tab=readme-ov-file). 

```@docs; canonical = false
SoilEnergyTemperatureClosure
```

## Methods

```@docs; canonical = false
get_thermal_properties
```

## Kernel functions

```@docs; canonical = false
compute_energy_tendency
```

```@docs; canonical = false
compute_energy_tendencies!
```

```@docs; canonical = false
compute_thermal_conductivity
```

## [References](@id "soilenergy.refs")

```@bibliography
Pages = ["soil_energy.md"]
Canonical = false
```
