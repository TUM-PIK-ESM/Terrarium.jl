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

variables(rre::RichardsEq) = (
    prognostic(:saturation_water_ice, XYZ(); closure = rre.closure, domain = UnitInterval(), desc = "Saturation level of water and ice in the pore space"),
    prognostic(:surface_excess_water, XY(), units = u"m", desc = "Excess water at the soil surface in m³/m²"),
    auxiliary(:water_table, XY(), units = u"m", desc = "Elevation of the water table in meters"),
    auxiliary(:hydraulic_conductivity, XYZ(z = Face()), units = u"m/s", desc = "Hydraulic conductivity of soil volumes in m/s"),
)

function initialize!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    closure!(state, model, get_closure(hydrology))
    return nothing
end

function compute_auxiliary!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    grid = get_grid(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    launch!(grid, :xyz, compute_hydraulic_conductivity!, state, grid, hydrology, strat, bgc)
    return nothing
end

function compute_tendencies!(state, model, hydrology::SoilHydrology{NF, <:RichardsEq}) where {NF}
    grid = get_grid(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    constants = get_constants(model)
    launch!(grid, :xyz, compute_saturation_tendency!, state, grid, hydrology, strat, bgc, constants)
    return nothing
end

# Kernels

# TODO: This is a dirty hack and basically the physical equivalent of money laundering.
# We should ideally use an implicit timestepping scheme to avoid this.
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
    @inbounds for k in 1:(N - 1)
        # TODO: This function might perform badly on GPU....
        # Can we optimize it somehow?
        if sat[i, j, k] > one(NF)
            # calculate excess saturation
            excess_sat = sat[i, j, k] - one(NF)
            # subtract excess water and add to layer above
            sat[i, j, k] -= excess_sat
            sat[i, j, k + 1] += excess_sat
        end
    end
    @inbounds if sat[i, j, N] > one(NF)
        # If the uppermost (surface) layer is oversaturated, add to excess water pool
        excess_sat = sat[i, j, N] - one(NF)
        sat[i, j, N] -= excess_sat
        state.surface_excess_water[i, j, 1] += excess_sat * Δzᵃᵃᶜ(i, j, N, field_grid)
    end
end

"""
    compute_water_table!(
        state,
        grid,
        ::SoilHydrology{NF},
        z_faces
    ) where {NF}

Kernel for diagnosing the water table at each grid point given the current soil saturation state. The argument
`z_faces` should be the z-coordinates of the grid on the layer faces.
"""
@kernel function compute_water_table!(
        state,
        grid,
        ::SoilHydrology{NF},
        z_faces
    ) where {NF}
    i, j = @index(Global, NTuple)
    sat = state.saturation_water_ice
    # scan z axis starting from the bottom (index 1) to find first non-saturated grid cell
    state.water_table[i, j, 1] = findfirst_z((i, j), <(one(NF)), z_faces, sat)
end

"""
    compute_hydraulic_conductivity!(
        state,
        grid,
        hydrology::SoilHydrology,
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry,
    )

Kernel for computing the hydraulic conductivity in all grid cells and soil layers.
"""
@kernel function compute_hydraulic_conductivity!(
        state,
        grid,
        hydrology::SoilHydrology,
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry,
    )
    i, j, k = @index(Global, NTuple)
    fgrid = get_field_grid(grid)
    if k <= 1
        state.hydraulic_conductivity[i, j, k] = hydraulic_conductivity(i, j, 1, fgrid, state, strat, hydrology, bgc)
    elseif k >= fgrid.Nz
        state.hydraulic_conductivity[i, j, k] = hydraulic_conductivity(i, j, fgrid.Nz, fgrid, state, strat, hydrology, bgc)
        state.hydraulic_conductivity[i, j, k + 1] = state.hydraulic_conductivity[i, j, k]
    else
        state.hydraulic_conductivity[i, j, k] = min_zᵃᵃᶠ(i, j, k, fgrid, hydraulic_conductivity, state, strat, hydrology, bgc)
    end
end

"""
    compute_saturation_tendency!(
        state,
        grid,
        hydrology::SoilHydrology,
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry,
    )

Kernel for computing the tendency of the prognostic `saturation_water_ice` variable in all grid cells and soil layers.
"""
@kernel function compute_saturation_tendency!(
        state,
        grid,
        hydrology::SoilHydrology,
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry,
        constants::PhysicalConstants
    )
    i, j, k = @index(Global, NTuple)
    # Get porosity
    por = porosity(i, j, k, state, hydrology, strat, bgc)
    # Compute volumetic water content tendency
    ∂θ∂t = volumetric_water_content_tendency(i, j, k, grid, state, hydrology, constants)
    # Rescale by porosity to get saturation tendency
    state.tendencies.saturation_water_ice[i, j, k] += ∂θ∂t / por
end

# Kernel functions

@inline function hydraulic_conductivity(i, j, k, grid, state, strat, hydrology, bgc)
    soil = soil_composition(i, j, k, state, strat, hydrology, bgc)
    return hydraulic_conductivity(hydrology.hydraulic_properties, soil)
end

@inline function volumetric_water_content_tendency(
        i, j, k, grid, state,
        hydrology::SoilHydrology{NF, <:RichardsEq},
        constants::PhysicalConstants
    ) where {NF}
    # Operators require the underlying Oceananigans grid
    field_grid = get_field_grid(grid)

    # Compute divergence of water fluxes
    # ∂θ∂t = ∇⋅K(θ)∇Ψ where Ψ = ψₘ + ψₕ + ψz + forcing (sources and sinks such as ET losses)
    ∂θ∂t = (
        - ∂zᵃᵃᶜ(i, j, k, field_grid, darcy_flux, state.pressure_head, state.hydraulic_conductivity)
            + forcing_ET(i, j, k, field_grid, state, hydrology.evapotranspiration, constants)
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
    Kₖ = (∇ψ < 0) * min(K[i, j, k - 1], K[i, j, k]) +
        (∇ψ >= 0) * min(K[i, j, k], K[i, j, k + 1])
    ## Note that elevation and hydrostatic pressure are assumed to be accounted
    ## for in the computation of ψ, so we do need any additional terms.
    q = -Kₖ * ∇ψ
    return q
end

# Matric potential <--> saturation closure relation

@kwdef struct SaturationPressureClosure <: AbstractClosureRelation end

closurevar(::SaturationPressureClosure) = auxiliary(:pressure_head, XYZ(), units = u"m", desc = "Total hydraulic pressure head in m water displaced at standard pressure")

function closure!(state, model::AbstractSoilModel, ::SaturationPressureClosure)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    z_faces = znodes(get_field_grid(grid), Center(), Center(), Face())
    z_centers = znodes(get_field_grid(grid), Center(), Center(), Center())
    # apply saturation correction
    launch!(grid, :xy, adjust_saturation_profile!, state, grid, hydrology)
    # update water table
    launch!(grid, :xy, compute_water_table!, state, grid, hydrology, z_faces)
    # determine pressure head from saturation
    launch!(grid, :xyz, saturation_to_pressure!, state, hydrology, strat, bgc, z_centers)
    return nothing
end

function invclosure!(state, model::AbstractSoilModel, ::SaturationPressureClosure)
    grid = get_grid(model)
    hydrology = get_soil_hydrology(model)
    strat = get_stratigraphy(model)
    bgc = get_biogeochemistry(model)
    z_faces = znodes(get_field_grid(grid), Center(), Center(), Face())
    z_centers = znodes(get_field_grid(grid), Center(), Center(), Center())
    # determine saturation from pressure
    launch!(grid, :xyz, pressure_to_saturation!, state, hydrology, strat, bgc, z_centers)
    # apply saturation correction
    launch!(grid, :xy, adjust_saturation_profile!, state, grid, hydrology)
    # update water table
    launch!(grid, :xy, compute_water_table!, state, grid, hydrology, z_faces)
    return nothing
end

@kernel function pressure_to_saturation!(
        state,
        hydrology::SoilHydrology{NF, <:RichardsEq},
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry,
        zs
    ) where {NF}
    i, j, k = @index(Global, NTuple)
    pressure_to_saturation!(i, j, k, state, hydrology, strat, bgc, zs)
end

@inline function pressure_to_saturation!(
        i, j, k, state,
        hydrology::SoilHydrology{NF, <:RichardsEq},
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry,
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
    por = porosity(i, j, k, state, hydrology, strat, bgc)
    vol_water_ice_content = swrc(ψm; θsat = por)
    state.saturation_water_ice[i, j, k] = vol_water_ice_content / por
    return nothing
end

@kernel function saturation_to_pressure!(
        state,
        hydrology::SoilHydrology{NF, <:RichardsEq},
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry,
        zs
    ) where {NF}
    i, j, k = @index(Global, NTuple)
    saturation_to_pressure!(i, j, k, state, hydrology, strat, bgc, zs)
end

@inline function saturation_to_pressure!(
        i, j, k, state,
        hydrology::SoilHydrology{NF, <:RichardsEq},
        strat::AbstractStratigraphy,
        bgc::AbstractSoilBiogeochemistry,
        zs
    ) where {NF}
    sat = state.saturation_water_ice[i, j, k] # assumed given
    # get inverse of SWRC
    inv_swrc = inv(get_swrc(hydrology))
    por = porosity(i, j, k, state, hydrology, strat, bgc)
    # compute matric pressure head
    ψm = inv_swrc(sat * por; θsat = por)
    # compute elevation pressure head
    ψz = zs[k]
    # compute hydrostatic pressure head assuming impermeable lower boundary
    z₀ = state.water_table[i, j, 1]
    ψh = max(0, z₀ - zs[k])
    # compute total pressure head as sum of ψh + ψm + ψz
    # note that ψh and ψz will cancel out in the saturated zone
    state.pressure_head[i, j, k] = ψh + ψm + ψz
    return nothing
end
