# Snow energy balance

```@meta
CurrentModule = Terrarium
```

```@setup snowenergy
using Terrarium
using InteractiveUtils
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

The prognostic snow energy is the depth-integrated (column) internal energy ``\bar{U}_\text{snow}`` (J/m²). Its tendency is computed as,
```math
\begin{equation}
\frac{\partial \bar{U}_\text{snow}}{\partial t} = Q_\text{base} - Q_\text{top} + Q_\text{precip} + Q_\text{subl},
\end{equation}
```
where ``Q_\text{top}`` is the net surface heat flux (the surface energy balance closure flux over the
snow), ``Q_\text{base}`` is the conductive heat flux at the snow base, ``Q_\text{precip}`` is the nondirectional (positive) sensible/latent heat advected by precipitation (see [Snow mass balance](@ref)), and ``Q_\text{subl}`` is an advective correction for energy loss due to sublimation. Meltwater drains as liquid water at 0 °C, which is the zero-enthalpy reference of the enthalpy closure below, so it carries no enthalpy out of the snowpack and thus no explicit meltwater energy flux is needed.

The sublimation correction arises from the same enthalpy reference. ``Q_\text{top}`` includes the surface latent heat flux which reduces the energy content of the snowpack by ``\rho_w L_{sg} E_\text{subl}``, but the the mass leaving the snowpack departs as ice, whose specific enthalpy relative to liquid water at 0 °C is ``-L_{sl}``. Adding back ``Q_\text{subl} = \rho_w L_{sl} E_\text{subl}`` leaves a net pack loss of ``\rho_w L_{lv} E_\text{subl}`` since ``L_{sg} = L_{sl} + L_{lv}`` by definition.

## Energy–temperature closure

The bulk snowpack is treated as an ice-water-air mixture, and its temperature and liquid water fraction are recovered from the depth-integrated internal energy using the same medium-agnostic [`FreeWater`](@extref FreezeCurves.FreeWater) enthalpy relations used for the soil (see [Soil energy balance](@ref)). The volumetric energy ``U_\text{snow} = \bar{U}_\text{snow} / d_\text{snow}`` is computed from the snow depth ``d_\text{snow}`` (see [`compute_volumetric_snow_energy`](@ref)). The `FreeWater` map references the internal energy to liquid water at 0 °C (``U_\text{snow} = 0``), so the phase-change range is ``U_\text{snow} \in [-L\theta, 0]`` with 0 °C ice at ``-L\theta`` where ``\theta`` is here the total volumetric content of water/ice. The volumetric latent heat is ``L\theta = \rho_\text{snow} L_{sl}`` — equivalently ``\rho_w L_{sl}`` times the water-equivalent content ``\theta = \rho_\text{snow}/\rho_w`` — so the density factor cancels in the depth integral, ``L\theta\, d_\text{snow} = \rho_w W_\text{snow} L_{sl}``, leaving the total latent heat set by the snow water equivalent alone. The bulk heat capacity ``C_\text{snow} = c_{p,i}\rho_i\theta_\text{ice} + c_{p,w}\rho_w\theta_\text{liq}`` is the ice/liquid volume-weighted sum (see [`compute_snow_volumetric_heat_capacity`](@ref)); the pore air enters only through the bulk density ``\rho_\text{snow}`` and its sensible heat is neglected. Because snow temperature cannot exceed 0 °C, the recovered free-water temperature is clipped at zero.

```@docs; canonical = false
SnowEnergyTemperatureClosure

compute_volumetric_snow_energy

compute_snow_volumetric_heat_capacity
```


## Advected heat from precipitation

Precipitation carries sensible and latent heat into the snowpack, referenced to liquid water at 0 °C (the
enthalpy reference of the closure, see [Snow energy balance](@ref)).

```@docs; canonical = false
compute_snow_precip_heat_flux
```

## Tendency

```@docs; canonical = false
compute_snow_energy_tendency
```

## Boundary and coupling heat fluxes

The snow→soil basal conductive flux and the blended soil-top flux, together with the snow surface
sublimation rate, are diagnosed after the surface energy balance solve.

```@docs; canonical = false
compute_snow_basal_heat_flux

compute_snow_soil_heat_flux

compute_snow_sublimation_flux

compute_snow_soil_boundary_fluxes!
```
