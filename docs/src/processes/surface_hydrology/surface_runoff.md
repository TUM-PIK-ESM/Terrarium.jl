# Surface runoff

```@meta
CurrentModule = Terrarium
```

```@setup runoff
using Terrarium
using InteractiveUtils
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

Surface runoff occurs when the rate of water reaching the surface exceeds the rate of infiltration into the soil. Surface runoff (overland flow) is often distinguished from other sources of subsurface runoff such as interflow in the vadose zone and baseflow in the saturated zone. Surface runoff is a key output of both hydrological and global land surface models since it serves as a input flux to lakes and rivers. Runoff that is routed into rivers can then serve as a key input to river routing schemes that produce freshwater fluxes for ocean components of Earth system models.

```@docs; canonical = false
AbstractSurfaceRunoff
```

```@example runoff
subtypes(Terrarium.AbstractSurfaceRunoff)
```

## Direct runoff

The simplest runoff scheme provided by Terrarium is [`DirectSurfaceRunoff`](@ref) which first routes all rainwater into the ground as infiltration and then simply routes the residual water flux as runoff.

Water in excess of instantaneous infiltration capacity accumulates as surface water and slowly drains back into the soil,
```math
\begin{equation}
D = \frac{S}{\tau_r}
\end{equation}
```
where $S$ is the accumulated surface water depth (m) and $\tau_r$ is the surface water drainage timescale (s).

```@docs; canonical = false
DirectSurfaceRunoff
```

```@example runoff
variables(DirectSurfaceRunoff(Float32))
```

### Infiltration capacity and soil saturation

Infiltration is limited by both the type and hydrological state of the soil. The effective infiltration can be calculated as
```math
\begin{equation}
I = \min(P_{\text{ground}}, I_{\max}) \times (1 - f_{\text{sat}})
\end{equation}
```
where $I_{\max}$ is here the maximum infiltration capacity determined by soil hydraulic conductivity (m/s) and $f_{\text{sat}}$ is the saturation fraction of the upper soil layer (0 = dry, 1 = saturated) (-).

When soil becomes saturated ($f_{\text{sat}} = 1$), infiltration drops to zero and all precipitation is routed to either surface runoff or surface water storage.

### Water table and drainage processes

This simplified surface runoff scheme does not explicitly model subsurface drainage or return flows from the water table. It is appropriate for upland areas where:
- Water table is deep relative to the surface
- Saturation excess runoff is minimal
- Hortonian overland flow (infiltration excess) is dominant

In riparian areas or regions with high water tables, more sophisticated approaches (e.g., variable source areas) would be necessary to capture saturation-excess runoff.

## Kernel functions

```@docs; canonical = false
compute_surface_drainage
```

```@docs; canonical = false
compute_infiltration
```

```@docs; canonical = false
compute_surface_runoff
```
