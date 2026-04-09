# Soil processes

```@meta
CurrentModule = Terrarium
```

```@setup coupled_soil
using Terrarium
```

!!! warning
    This page is a work in progress. If you have any questions or notice any errors, please [raise an issue](https://github.com/NumericalEarth/Terrarium.jl/issues).

Subtypes of [`AbstractSoil`](@ref) represent a set of coupled processes affecting natural soils. The three main soil processes currently considered in Terrarium are:

- **Energy** — heat conduction with phase change of water/ice in soil pores
- **Hydrology** — hydrological properties and dynamics of the soil
- **Biogeochemistry** — soil carbon cycling and biogeochemical processes

```@docs; canonical = false
SoilEnergyWaterCarbon
```

```@example coupled_soil
variables(SoilEnergyWaterCarbon(Float32))
```

## Process Interface

```@docs; canonical = false
initialize!(state, grid, ::SoilEnergyWaterCarbon, constants::PhysicalConstants)
compute_auxiliary!(state, grid, ::SoilEnergyWaterCarbon, constants::PhysicalConstants)
compute_tendencies!(state, grid, ::SoilEnergyWaterCarbon, constants::PhysicalConstants)
```

## Methods

```@docs; canonical = false
get_stratigraphy(soil::AbstractSoil)
get_energy_balance(soil::AbstractSoil)
get_hydrology(soil::AbstractSoil) 
get_biogeochemistry(soil::AbstractSoil)
```

## Component processes

See the following pages for more details on the various implementations and relevant parameterizations for each soil process:

- [Soil energy balance](@ref) — heat conduction and freeze-thaw
- [Soil hydrology](@ref) — water transport and saturation
- [Soil stratigraphy](@ref) — vertical distribution of material properties