"""
    RichardsEq{PS} <: AbstractVerticalFlow

Operator for soil hydrology corresponding to the Richardson-Richards equation for variably saturated
flow in porous media.
"""
@kwdef struct RichardsEq{PS} <: AbstractVerticalFlow
    "Closure relation for mapping between water potential (hydraulic head) and saturation"
    closure::PS = SaturationPressureClosure()
end

get_closure(hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF} = hydrology.vertflow.closure

variables(hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF} = (
    prognostic(:saturation_water_ice, XYZ(); closure=get_closure(hydrology), domain=UnitInterval(), desc="Saturation level of water and ice in the pore space"),
    prognostic(:surface_excess_water, XY(), units=u"m", desc="Excess water at the soil surface in m³/m²"),
    auxiliary(:hydraulic_conductivity, XYZ(z=Face()), units=u"m/s", desc="Hydraulic conductivity of soil volumes in m/s"),
    auxiliary(:water_table, XY(), units=u"m", desc="Elevation of the water table in meters"),
    input(:liquid_water_fraction, XYZ(), default=NF(1), domain=UnitInterval(), desc="Fraction of unfrozen water in the pore space"), 
)

function initialize!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    set!(state.liquid_water_fraction, 1)
    closure!(state, model, get_closure(hydrology))
    return nothing
end

function compute_auxiliary!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    grid = get_grid(model)
    strat = get_soil_stratigraphy(model)
    bgc = get_soil_biogeochemistry(model)
    launch!(state, grid, :xyz, compute_hydraulics_kernel!, hydrology, strat, bgc)
    return nothing
end

function compute_tendencies!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    grid = get_grid(model)
    strat = get_soil_stratigraphy(model)
    constants = get_constants(model)
    launch!(state, grid, :xyz, compute_saturation_tendency!, hydrology, strat, constants, nothing)
    return nothing
end

function compute_tendencies!(state, model::AbstractLandModel, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    grid = get_grid(model)
    strat = get_soil_stratigraphy(model)
    evapotranspiration = get_evapotranspiration(model)
    constants = get_constants(model)
    launch!(state, grid, :xyz, compute_saturation_tendency!, hydrology, strat, constants, evapotranspiration)
    return nothing
end

# Kernels

"""
    adjust_saturation_profile!(
        state,
        grid,
        ::SoilHydrology{NF}
    )

Kernel for adjusting saturation profiles to account for oversaturation due to numerical error.
This implementation scans over the saturation profiles at each lateral grid cell and redistributes
excess water upward layer-by-layer until reaching the topmost layer, where any remaining excess
water is added to the `surface_excess_water` pool.
"""
@kernel function adjust_saturation_profile!(
    state,
    grid,
    ::SoilHydrology{NF}
) where {NF}
    i, j = @index(Global, NTuple)
    sat = state.saturation_water_ice
    field_grid = get_field_grid(grid)
    N = field_grid.Nz
    # First iterate over soil layers from bottom to top
    # TODO: This function might perform badly on GPU....
    # Can we optimize it somehow?
    @inbounds for k in 1:N-1
        if sat[i, j, k] > one(NF)
            # calculate excess saturation
            excess_sat = sat[i, j, k] - one(NF)
            # subtract excess water and add to layer above;
            # note that we need to rescale by the cell thickness to properly conserve mass
            sat[i, j, k] -= excess_sat
            sat[i, j, k+1] += excess_sat * Δzᵃᵃᶜ(i, j, k, field_grid) / Δzᵃᵃᶜ(i, j, k+1, field_grid)
        end
    end
    # then from top to bottom
    @inbounds for k in N:-1:2
        if sat[i, j, k] < zero(NF)
            # calculate saturation deficit
            deficit_sat = -sat[i, j, k]
            # add back saturation deficit and subtract from layer below
            sat[i, j, k] += deficit_sat
            sat[i, j, k-1] -= deficit_sat * Δzᵃᵃᶜ(i, j, k, field_grid) / Δzᵃᵃᶜ(i, j, k-1, field_grid)
        end
    end
    @inbounds if sat[i, j, N] > one(NF)
        # If the uppermost (surface) layer is oversaturated, add to excess water pool
        excess_sat = sat[i, j, N] - one(NF)
        sat[i, j, N] -= excess_sat
        state.surface_excess_water[i, j, 1] += excess_sat * Δzᵃᵃᶜ(i, j, N, field_grid)
    end
    @inbounds if sat[i, j, 1] < zero(NF)
        # If the uppermost (surface) layer has a deficit, just set to zero.
        # This constitutes a mass balance violation but should not happen under realistic conditions.
        sat[i, j, 1] = zero(NF)
    end
end

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
    state,
    grid,
    hydrology::SoilHydrology{NF, <:RichardsEq},
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
    compute_saturation_tendency!(
        state,
        grid,
        hydrology::SoilHydrology,
        strat::AbstractStratigraphy,
    )

Kernel for computing the tendency of the prognostic `saturation_water_ice` variable in all grid cells and soil layers.
"""
@kernel function compute_saturation_tendency!(
    state,
    grid,
    hydrology::SoilHydrology,
    strat::AbstractStratigraphy,
    constants::PhysicalConstants,
    evapotranspiration::Union{Nothing, AbstractEvapotranspiration}
)
    i, j, k = @index(Global, NTuple)
    # Compute volumetic water content tendency
    ∂θ∂t = volumetric_water_content_tendency(i, j, k, grid, state, hydrology, constants, evapotranspiration)
    # Get porosity
    por = porosity(i, j, k, state, grid, strat)
    # Rescale by porosity to get saturation tendency
    state.tendencies.saturation_water_ice[i, j, k] +=  ∂θ∂t / por
end

# Kernel functions

# This function is needed for an Oceananigans grid operator
@inline function hydraulic_conductivity(i, j, k, grid, state, hydrology, strat, bgc)
    soil = soil_volume(i, j, k, state, grid, strat, hydrology, bgc)
    return hydraulic_conductivity(hydrology.hydraulic_properties, soil)
end

@inline function volumetric_water_content_tendency(
    i, j, k, grid, state,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    constants::PhysicalConstants,
    evapotranspiration::Union{Nothing, AbstractEvapotranspiration}
) where {NF}
    # Operators require the underlying Oceananigans grid
    field_grid = get_field_grid(grid)

    # Compute divergence of water fluxes
    # ∂θ∂t = ∇⋅K(θ)∇Ψ + forcing, where Ψ = ψₘ + ψₕ + ψz, and "forcing" represents sources and sinks such as ET losses
    ∂θ∂t = (
        - ∂zᵃᵃᶜ(i, j, k, field_grid, darcy_flux, state.pressure_head, state.hydraulic_conductivity)
        + forcing(i, j, k, state, grid, evapotranspiration, hydrology, constants)
        + forcing(i, j, k, state, grid, hydrology.vwc_forcing, hydrology)
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

# Matric potential <--> saturation closure relation

@kwdef struct SaturationPressureClosure <: AbstractClosureRelation end

closurevar(::SaturationPressureClosure) = auxiliary(:pressure_head, XYZ(), units=u"m", desc="Total hydraulic pressure head in m water displaced at standard pressure")

function closure!(state, model::AbstractSoilModel, ::SaturationPressureClosure)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_soil_stratigraphy(model)
    z_centers = znodes(get_field_grid(grid), Center(), Center(), Center())
    # apply saturation correction
    launch!(state, grid, :xy, adjust_saturation_profile!, hydrology)
    # update water table
    compute_water_table!(state, grid, hydrology)
    # determine pressure head from saturation
    launch!(state, grid, :xyz, saturation_to_pressure!, hydrology, strat, z_centers)
    return nothing
end

function invclosure!(state, model::AbstractSoilModel, ::SaturationPressureClosure)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_soil_stratigraphy(model)
    z_centers = znodes(get_field_grid(grid), Center(), Center(), Center())
    # determine saturation from pressure
    launch!(state, grid, :xyz, pressure_to_saturation!, hydrology, strat, z_centers)
    # apply saturation correction
    launch!(state, grid, :xy, adjust_saturation_profile!, hydrology)
    # update water table
    compute_water_table!(state, grid, hydrology)
    return nothing
end

@kernel function pressure_to_saturation!(
    state, grid,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    strat::AbstractStratigraphy,
    zs
) where {NF}
    i, j, k = @index(Global, NTuple)
    pressure_to_saturation!(i, j, k, state, grid, hydrology, strat, zs)
end

@kernel function saturation_to_pressure!(
    state, grid,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    strat::AbstractStratigraphy,
    zs
) where {NF}
    i, j, k = @index(Global, NTuple)
    saturation_to_pressure!(i, j, k, state, grid, hydrology, strat, zs)
end

@inline function pressure_to_saturation!(
    i, j, k, state, grid,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    strat::AbstractStratigraphy,
    zs
) where {NF}
    ψ = state.pressure_head[i, j, k] # assumed given
    # compute elevation pressure head
    ψz = zs[k]
    # compute hydrostatic pressure head assuming impermeable lower boundary
    z₀ = state.water_table[i, j, 1]
    ψh = max(0, z₀ - zs[k])
    # remove hydrostatic and elevation components
    ψm = ψ - ψh - ψz
    swrc = get_swrc(hydrology)
    por = porosity(i, j, k, state, grid, strat)
    vol_water_ice_content = swrc(ψm; θsat=por)
    state.saturation_water_ice[i, j, k] = vol_water_ice_content / por
    return nothing
end

@inline function saturation_to_pressure!(
    i, j, k, state, grid,
    hydrology::SoilHydrology{NF, <:RichardsEq},
    strat::AbstractStratigraphy,
    zs
) where {NF}
    sat = state.saturation_water_ice[i, j, k] # assumed given
    # get inverse of SWRC
    inv_swrc = inv(get_swrc(hydrology))
    por = porosity(i, j, k, state, grid, strat)
    # compute matric pressure head
    ψm = inv_swrc(sat*por; θsat=por)
    # compute elevation pressure head
    ψz = zs[k]
    # compute hydrostatic pressure head assuming impermeable lower boundary
    # TODO: can we generalize this for arbitrary lower boundaries?
    z₀ = state.water_table[i, j, 1]
    ψh = max(0, z₀ - zs[k])
    # compute total pressure head as sum of ψh + ψm + ψz
    # note that ψh and ψz will cancel out in the saturated zone
    state.pressure_head[i, j, k] = ψh + ψm + ψz
    return nothing
end
