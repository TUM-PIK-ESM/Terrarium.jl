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

The prognostic snow energy is the depth-integrated (column) internal energy ``U_\text{snow}`` (J/m²). Its tendency is computed as,
```math
\begin{equation}
\frac{\partial U_\text{snow}}{\partial t} = Q_\text{base} - Q_\text{top} + Q_\text{precip} + Q_\text{subl},
\end{equation}
```
where ``Q_\text{top}`` is the net surface heat flux (the surface energy balance closure flux over the
snow), ``Q_\text{base}`` is the conductive heat flux at the snow base, ``Q_\text{precip}`` is the
sensible/latent heat advected by precipitation (see [Snow mass balance](@ref)), and ``Q_\text{subl}`` is
an advective correction for energy loss due to sublimation. Meltwater drains as liquid water at 0 °C, which is the zero-enthalpy reference of the enthalpy closure below, so it carries no enthalpy out of the snowpack and no
explicit meltwater energy term appears.

The sublimation correction arises from the same enthalpy reference.``Q_\text{top}`` includes the surface latent heat flux which reduces the energy content of the snowpack by ``\rho_w L_s E_\text{subl}``, but the the mass leaving the snowpack departs as ice, whose specific enthalpy relative to liquid water at 0 °C is ``-L_f``. Adding back ``Q_\text{subl} = \rho_w L_f E_\text{subl}`` leaves a net pack loss of ``\rho_w L_v E_\text{subl}`` since ``L_s = L_f + L_v`` by definition.

## Energy–temperature closure

The bulk snowpack is treated as an ice-water-air mixture, and its temperature and liquid water fraction are recovered from the depth-integrated internal energy using the same medium-agnostic [`FreeWater`](@extref FreezeCurves.FreeWater) enthalpy relations used for the soil (see [Soil energy balance](@ref)). The volumetric energy ``U_v = E / d_s`` is computed from the snow depth ``d_s`` (see [`volumetric_snow_energy`](@ref)) without allocating an auxiliary field. The `FreeWater` map references the internal energy to liquid water at 0 °C (``U = 0``), so the phase-change band is ``U_v \in [-L\theta, 0]`` with ice at 0 °C at ``-L\theta``. Because snow temperature cannot exceed 0 °C, the recovered free-water temperature is clipped at zero.

```@docs; canonical = false
SnowEnergyTemperatureClosure

volumetric_snow_energy
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
