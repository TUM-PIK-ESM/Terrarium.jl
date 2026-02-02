"""
    RichardsEq{PS} <: AbstractVerticalFlow

[`SoilHydrology`](@ref) flow operator implementing the mixed saturation-pressure form
of the Richardson-Richards equation:

```math
\\phi(z) \\frac{\\partial s_{\\mathrm{wi}}(\\Psi(z,t))}{\\partial t} = \\n+\\nabla \\cdot \\bigl(K(s_{\\mathrm{wi}}, T) \\; \\n+\\nabla (\\psi_m + 1)\\bigr)
```
which describes the vertical movement of water according to gravity-driven
percolation and capillary-driven diffusion.

State variables defined by the Richards' formulation of `SoilHydrology`:

- `saturation_water_ice`: saturation level of water and ice in the pore space.
- `surface_excess_water`: excess water at the soil surface (m^3/m^2).
- `hydraulic_conductivity`: hydraulic conductivity at cell centers (m/s).
- `water_table``: elevation of the water table (m).
- `liquid_water_fraction`: fraction of unfrozen liquid water in the pore space (dimensionless).

See also [`SaturationPressureClosure`](@ref) and [`SoilHydraulics`](@ref) for details regarding the
closure relating saturtion and pressure head.
"""
@kwdef struct RichardsEq <: AbstractVerticalFlow end

variables(hydrology::SoilHydrology{NF, RichardsEq}) where {NF} = (
    prognostic(:saturation_water_ice, XYZ(); closure=get_closure(hydrology), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
    prognostic(:surface_excess_water, XY(), units=u"m", desc="Excess water at the soil surface in m³/m²"),
    auxiliary(:hydraulic_conductivity, XYZ(z=Face()), units=u"m/s", desc="Hydraulic conductivity of soil volumes in m/s"),
    auxiliary(:water_table, XY(), units=u"m", desc="Elevation of the water table in meters"),
    input(:liquid_water_fraction, XYZ(), default=NF(1), domain=UnitInterval(), desc="Fraction of unfrozen water in the pore space"), 
)

@propagate_inbounds surface_excess_water(i, j, grid, state, ::SoilHydrology{NF, RichardsEq}) where {NF} = state.surface_excess_water[i, j]

function initialize!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    constants::PhysicalConstants,
    args...
) where {NF}
    closure!(state, grid, hydrology, soil)
    compute_auxiliary!(state, grid, hydrology, soil, constants)
    return nothing
end

function compute_auxiliary!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    args...
) where {NF}
    strat = get_stratigraphy(soil)
    bgc = get_biogeochemistry(soil)
    out = auxiliary_fields(state, hydrology)
    fields = get_fields(state, hydrology, bgc; except=out)
    launch!(grid, XYZ, compute_hydraulics_kernel!, out, fields, hydrology, strat, bgc)
    return nothing
end

function compute_tendencies!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    constants::PhysicalConstants,
    evtr::Optional{AbstractEvapotranspiration} = nothing,
    args...
) where {NF}
    strat = get_stratigraphy(soil)
    bgc = get_biogeochemistry(soil)
    tendencies = tendency_fields(state, hydrology)
    fields = get_fields(state, hydrology, bgc, evtr)
    clock = state.clock
    launch!(grid, XYZ, compute_saturation_tendency_kernel!,
            tendencies, clock, fields, hydrology, strat, bgc, constants, evtr)
    return nothing
end

"""
    $TYPEDSIGNATURES

Computes `pressure_head` ``Ψ = ψm + ψz + ψh`` from the current `saturation_water_ice` state.
"""
function closure!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    args...
) where {NF}
    # apply saturation correction
    adjust_saturation_profile!(state, grid, hydrology)
    # update water table
    compute_water_table!(state, grid, hydrology)
    # determine pressure head from saturation
    strat = get_stratigraphy(soil)
    bgc = get_biogeochemistry(soil)
    out = (pressure_head = state.pressure_head,)
    fields = get_fields(state, hydrology, bgc; except = out)
    launch!(grid, XYZ, saturation_to_pressure_kernel!,
            out, fields, hydrology.closure, hydrology, strat, bgc)
    return nothing
end

"""
    $TYPEDSIGNATURES

Computes `saturation_water_ice` from the current `pressure_head` state.
"""
function invclosure!(
    state, grid,
    hydrology::SoilHydrology{NF, RichardsEq},
    soil::AbstractSoil,
    args...
) where {NF}
    strat = get_stratigraphy(soil)
    bgc = get_biogeochemistry(soil)
    out = (saturation_water_ice = state.saturation_water_ice,)
    fields = get_fields(state, hydrology, bgc; except = out)
    # determine saturation from pressure
    launch!(grid, XYZ, pressure_to_saturation_kernel!,
            out, fields, hydrology.closure, hydrology, strat, bgc)
    # apply saturation correctionh
    adjust_saturation_profile!(state, grid, hydrology)
    # update water table
    compute_water_table!(state, grid, hydrology)
    return nothing
end


# Kernels

"""
    $TYPEDSIGNATURES

Kernel function for computing the tendency of the prognostic `saturation_water_ice` variable in all grid cells and soil layers.
"""
@propagate_inbounds function compute_saturation_tendency!(
    saturation_water_ice_tendency, i, j, k, grid, clock, fields,
    hydrology::SoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
    evapotranspiration::Union{Nothing, AbstractEvapotranspiration}
)
    # Compute volumetic water content tendency
    ∂θ∂t = volumetric_water_content_tendency(i, j, k, grid, clock, fields, hydrology, constants, evapotranspiration)
    # Get porosity
    por = porosity(i, j, k, grid, fields, strat, bgc)
    # Rescale by porosity to get saturation tendency
    saturation_water_ice_tendency[i, j, k] +=  ∂θ∂t / por
end

# Kernel functions

"""
    $SIGNATURES

Compute the volumetric water content (VWC) tendency at grid cell `i, j k`. Note that the
VWC tendency is not scaled by the porosity and is thus not a saturation tendency.
"""
@inline function volumetric_water_content_tendency(
    i, j, k, grid, clock, fields,
    hydrology::SoilHydrology{NF, RichardsEq},
    constants::PhysicalConstants,
    evapotranspiration::Union{Nothing, AbstractEvapotranspiration}
) where {NF}
    # Operators require the underlying Oceananigans grid
    field_grid = get_field_grid(grid)

    # Compute divergence of water fluxes
    # ∂θ∂t = ∇⋅K(θ)∇Ψ + forcing, where Ψ = ψₘ + ψₕ + ψz, and "forcing" represents sources and sinks such as ET losses
    ∂θ∂t = (
        - ∂zᵃᵃᶜ(i, j, k, field_grid, darcy_flux, fields.pressure_head, fields.hydraulic_conductivity)
        + forcing(i, j, k, grid, clock, fields, evapotranspiration, hydrology, constants) # ET forcing
        + forcing(i, j, k, grid, clock, fields, hydrology.vwc_forcing, hydrology) # generic user-defined forcing
    )
    return ∂θ∂t
end

"""
    $SIGNATURES

Kernel function for computing the Darcy flux over layer faces from the pressure head `ψ` and hydraulic
conductivity `K`.
"""
@inline function darcy_flux(i, j, k, grid, ψ, K)
    # Darcy's law: q = -K ∂ψ/∂z
    # TODO: also account for effect of temperature on matric potential
    ∇ψ = ∂zᵃᵃᶠ(i, j, k, grid, ψ)
    ## Take minimum of hydraulic conductivities in the direction of flow
    Kₖ = (∇ψ < 0)*min(K[i, j, k-1], K[i, j, k]) +
         (∇ψ >= 0)*min(K[i, j, k], K[i, j, k+1])
    ## Note that elevation and hydrostatic pressure are assumed to be accounted
    ## for in the computation of ψ, so we don't need any additional terms.
    q = -Kₖ*∇ψ
    return q
end

"""
    $SIGNATURES

Compute the hydraulic conductivity at the center of the grid cell `i, j, k`.
"""
@inline function hydraulic_conductivity(
    i, j, k, grid, state,
    hydrology::SoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    soil = soil_volume(i, j, k, grid, state, strat, hydrology, bgc)
    return hydraulic_conductivity(hydrology.hydraulic_properties, soil)
end

# Kernels

@kernel inbounds=true function compute_saturation_tendency_kernel!(
    tendencies, grid, clock, fields,
    hydrology::SoilHydrology{NF, RichardsEq},
    args...
) where {NF}
    i, j, k = @index(Global, NTuple)
    compute_saturation_tendency!(tendencies.saturation_water_ice, i, j, k, grid, clock, fields, hydrology, args...)
end
