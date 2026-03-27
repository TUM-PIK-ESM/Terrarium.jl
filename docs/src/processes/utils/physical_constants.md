```@meta
CurrentModule = Terrarium
```

```@setup default
using Terrarium
```

# Physical constants

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

## Overview

`PhysicalConstants` collects fundamental physical constants used throughout Terrarium's
process implementations. All constants are stored as fields of a single struct so that
they are passed explicitly through the call graph — avoiding global state and keeping
the code fully differentiable with Enzyme.jl. The struct is parametrically typed so that
constants are automatically promoted to the model's numeric precision `NF`.

Default values follow standard references. Individual constants can be overridden at
construction to support unit-testing or sensitivity studies.

```@docs; canonical = false
PhysicalConstants
```

```@example default
PhysicalConstants(Float32)
```

## Methods

```@docs; canonical = false
celsius_to_kelvin
```

```@docs; canonical = false
stefan_boltzmann
```

```@docs; canonical = false
psychrometric_constant
```

```@docs; canonical = false
compute_vpd
```
