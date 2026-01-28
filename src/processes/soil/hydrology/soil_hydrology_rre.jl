"""
    RichardsEq{PS} <: AbstractVerticalFlow

[`SoilHydrology`](@ref) flow operator implementing the mixed saturation-pressure form
of the Richardson-Richards equation:

```math
\\phi(z) \\frac{\\partial s_{\\mathrm{wi}}(\\Psi(z,t))}{\\partial t} = \\n+\\nabla \\cdot \\bigl(K(s_{\\mathrm{wi}}, T) \\; \\n+\\nabla (\\psi_m + 1)\\bigr)
```
which describes the vertical movement of water according to gravity-driven
percolation and capillary-driven diffusion.

State variables defined by this operator:

- `saturation_water_ice`: saturation level of water and ice in the pore space.
- `surface_excess_water`: excess water at the soil surface (m^3/m^2).
- `hydraulic_conductivity`: hydraulic conductivity at cell centers (m/s).
- `water_table``: elevation of the water table (m).
- `liquid_water_fraction`: fraction of unfrozen liquid water in the pore space (dimensionless).

See also [`SaturationPressureClosure`](@ref) and [`SoilHydraulics`](@ref) for details regarding the
closure relating saturtion and pressure head.
"""
@kwdef struct RichardsEq <: AbstractVerticalFlow end

"""
Alias for `SoilHydrology{NF, RichardsEq, Closure}`
"""
const SoilHydrologyRRE{NF, Closure} = SoilHydrology{NF, RichardsEq, Closure}

variables(hydrology::SoilHydrologyRRE{NF}) where {NF} = (
    prognostic(:saturation_water_ice, XYZ(); closure=get_closure(hydrology), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
    prognostic(:surface_excess_water, XY(), units=u"m", desc="Excess water at the soil surface in m³/m²"),
    auxiliary(:hydraulic_conductivity, XYZ(z=Face()), units=u"m/s", desc="Hydraulic conductivity of soil volumes in m/s"),
    auxiliary(:water_table, XY(), units=u"m", desc="Elevation of the water table in meters"),
    input(:liquid_water_fraction, XYZ(), default=NF(1), domain=UnitInterval(), desc="Fraction of unfrozen water in the pore space"), 
)

@propagate_inbounds surface_excess_water(i, j, grid, state, ::SoilHydrologyRRE{NF}) where {NF} = state.surface_excess_water[i, j]

function initialize!(
    state, grid,
    hydrology::SoilHydrologyRRE{NF},
    strat::AbstractStratigraphy
) where {NF}
    closure!(state, grid, hydrology, strat)
    return nothing
end

function compute_auxiliary!(
    state, grid,
    hydrology::SoilHydrologyRRE{NF},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    launch!(grid, XYZ, compute_hydraulics_kernel!, state, hydrology, strat, bgc)
    return nothing
end

function compute_tendencies!(
    state, grid,
    hydrology::SoilHydrologyRRE{NF},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants
) where {NF}
    launch!(grid, XYZ, compute_saturation_tendency_kernel!, state, hydrology, strat, bgc, constants, nothing)
    return nothing
end

function compute_tendencies!(
    state, grid,
    hydrology::SoilHydrologyRRE{NF},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
    evtr::AbstractEvapotranspiration
) where {NF}
    launch!(grid, XYZ, compute_saturation_tendency_kernel!, state, hydrology, strat, bgc, constants, evtr)
    return nothing
end

# Kernels

"""
    compute_hydraulics_kernel!(
        state,
        grid,
        hydrology::SoilHydrology,
        strat::AbstractStratigraphy
    )

Kernel for computing soil hydraulics and unsaturated hydraulic conductivity.
"""
@kernel function compute_hydraulics_kernel!(
    state, grid,
    hydrology::SoilHydrologyRRE{NF},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    i, j, k = @index(Global, NTuple)
    fgrid = get_field_grid(grid)

    # compute hydraulic conductivity
    @inbounds if k <= 1
        state.hydraulic_conductivity[i, j, k] = hydraulic_conductivity(i, j, 1, fgrid, state, hydrology, strat, bgc)
    elseif k >= fgrid.Nz
        state.hydraulic_conductivity[i, j, k] = hydraulic_conductivity(i, j, fgrid.Nz, fgrid, state, hydrology, strat, bgc)
        state.hydraulic_conductivity[i, j, k + 1] = state.hydraulic_conductivity[i, j, k]
    else
        state.hydraulic_conductivity[i, j, k] = min_zᵃᵃᶠ(i, j, k, fgrid, hydraulic_conductivity, state, hydrology, strat, bgc)
    end
end

"""
    compute_saturation_tendency_kernel!(
        state,
        grid,
        hydrology::SoilHydrology,
        strat::AbstractStratigraphy,
    )

Kernel for computing the tendency of the prognostic `saturation_water_ice` variable in all grid cells and soil layers.
"""
@kernel function compute_saturation_tendency_kernel!(
    state, grid,
    hydrology::SoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry,
    constants::PhysicalConstants,
    evapotranspiration::Union{Nothing, AbstractEvapotranspiration}
)
    i, j, k = @index(Global, NTuple)
    # Compute volumetic water content tendency
    ∂θ∂t = volumetric_water_content_tendency(i, j, k, grid, state, hydrology, constants, evapotranspiration)
    # Get porosity
    por = porosity(i, j, k, grid, state, strat, bgc)
    # Rescale by porosity to get saturation tendency
    state.tendencies.saturation_water_ice[i, j, k] +=  ∂θ∂t / por
end

# Kernel functions

"""
    $SIGNATURES

Compute the volumetric water content (VWC) tendency at grid cell `i, j k`. Note that the
VWC tendency is not scaled by the porosity and is thus not a saturation tendency.
"""
@inline function volumetric_water_content_tendency(
    i, j, k, grid, state,
    hydrology::SoilHydrologyRRE{NF},
    constants::PhysicalConstants,
    evapotranspiration::Union{Nothing, AbstractEvapotranspiration}
) where {NF}
    # Operators require the underlying Oceananigans grid
    field_grid = get_field_grid(grid)

    # Compute divergence of water fluxes
    # ∂θ∂t = ∇⋅K(θ)∇Ψ + forcing, where Ψ = ψₘ + ψₕ + ψz, and "forcing" represents sources and sinks such as ET losses
    ∂θ∂t = (
        - ∂zᵃᵃᶜ(i, j, k, field_grid, darcy_flux, state.pressure_head, state.hydraulic_conductivity)
        + forcing(i, j, k, grid, state, evapotranspiration, hydrology, constants)
        + forcing(i, j, k, grid, state, hydrology.vwc_forcing, hydrology)
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
    ## for in the computation of ψ, so we do need any additional terms.
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
