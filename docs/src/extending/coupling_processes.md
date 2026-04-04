# Coupling processes

```@meta
CurrentModule = Terrarium
```

```@setup coupling
using Terrarium
using Terrarium: auxiliary, input
```

Most realistic land surface simulations require multiple processes that must share information at each time step — for example, a soil energy balance that needs the surface skin temperature computed by a skin temperature scheme, or a canopy interception scheme that consumes the leaf area index produced by the vegetation module. Terrarium allows for two distinct coupling strategies that can be chosen depending on how tightly the coupled processes are bound to one another.

## Weak coupling

The simplest way to couple processes together in Terrarium is via [`input`](@ref) variables. This strategy takes advantage the variable promotion rules described on the [State variables](@ref) page.

1. Process A declares a variable `X` as an `input`. Input variables are agnostic to their source and should always be treated as strictly read-only.
2. Process B declares the same variable as `prognostic` or `auxiliary`.
3. When the variables of both processes are combined in a [`Variables`](@ref) container, the `input` declaration from A is automatically promoted to reference B's `prognostic`/`auxiliary` field. A single `Field` is allocated for the shared variable; both processes read from it without either one knowing about the other's type.

The advantage of weak coupling is that the implementations of processes A and B remain completely **independent** of each other: neither carries a type parameter for the other, and their kernel functions only reference `fields.X` as an ordinary `Field`. This means both processes can be developed, tested, and reused in isolation.

### Example: Ground surface temperature

A concrete example in Terrarium is the `ground_temperature` variable which represents the temperature of the uppermost subsurface layer. `ground_temperature` is defined as a derived auxiliary variable by [`SoilEnergyBalance`](@ref), with the resulting `Field` being simply a view of the uppermost soil `temperature` layer:

```julia
variables(energy::SoilEnergyBalance) = (
    prognostic(:internal_energy, XYZ(); closure = energy.closure, units = u"J/m^3", desc = "Internal energy of the soil volume, including both latent and sensible components"),
    auxiliary(:ground_temperature, XY(), ground_temperature, energy, units = u"°C", desc = "Temperature of the uppermost ground or soil grid cell in °C"),
)

function ground_temperature(energy::SoilEnergyBalance, grid, clock, fields)
    fgrid = get_field_grid(grid)
    # Use uppermost soil layer as ground temperature
    return @view fields.temperature[:, :, fgrid.Nz]
end
```

This `ground_temperature` is then consumed by several other processes such as the [Surface Energy Balance](@ref), [canopy evapotranspiration](@ref "Canopy evapotranspiration"), and the vegetation stress factors computed for [autotrophic respiration](@ref "Autotrophic respiration"). These processes can simply declare `ground_temperature` as an input variable:

```julia
input(:ground_temperature, XY(), default = 10.0, units = u"°C")
```

When used standalone, `ground_temperature` will be allocated as an independent input `Field` (here with a uniform initial value of 10°C across space). This allows for easy testing in isolation.

## Strong coupling

Strong coupling is appropriate when two or more processes are tightly intertwined; this may be due to the processes sharing kernel function dispatches or requiring access to each other's parameters.

TODO

## Coupled process types

```@docs; canonical = false
AbstractCoupledProcesses
```

The standard pattern is to define a `struct` whose fields are the sub-processes, each with a concrete abstract type bound as a type parameter:

```julia
@kwdef struct MyCoupledProcess{
        NF,
        ProcA <: AbstractProcessA{NF},
        ProcB <: AbstractProcessB{NF},
    } <: AbstractCoupledProcesses{NF}
    procA::ProcA = ConcreteProcessA(NF)
    procB::ProcB = ConcreteProcessB(NF)
end
```

Kernel functions for the coupled process can then dispatch on or delegate to the sub-process types:

```julia
@propagate_inbounds function compute_coupled_flux!(
        tendencies, i, j, k, grid, fields,
        proc::MyCoupledProcess,
        args...
    )
    # Access sub-process parameters directly
    α = proc.proc_a.some_parameter

    # Delegate to sub-process kernel functions
    flux_a = compute_flux_a(i, j, k, grid, fields, proc.proc_a)
    flux_b = compute_flux_b(i, j, k, grid, fields, proc.proc_b, α)
    tendencies.X[i, j, k] += flux_a + flux_b
    return tendencies
end
```

### Example: `SurfaceEnergyBalance`

[`SurfaceEnergyBalance`](@ref) is a concrete `AbstractCoupledProcesses` that composes four sub-processes, [`AbstractSkinTemperature`](@ref), [`AbstractTurbulentFluxes`](@ref), [`AbstractRadiativeFluxes`](@ref), and [`AbstractAlbedo`](@ref), into a single surface energy calculation.

Its joint kernel function passes each sub-process to the appropriate kernel function, dispatching on their concrete types to select the correct implementation.

This pattern allows every combination of sub-process implementations to be statically compiled into a specialized method at zero runtime overhead, while keeping each sub-process's physics fully encapsulated in its own kernel functions.
