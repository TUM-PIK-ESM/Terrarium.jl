# Saturation ↔ Hydraulic (pressure) head

"""
    $TYPEDEF

Represents a closure relating saturation of water/ice in soil pores to a corresponding
pressure (or hydraulic) head. Note that here "pressure head" is defined to be synonymous
with hydraulic head, i.e. including all both elevation and hydrostatic pressure contributions.
This relation is typically described by soil property-dependent *soil-water retention curve* (SWRC)
which is here defined in [`SoilHydraulics`](@ref).
"""
@kwdef struct SaturationPressureClosure <: AbstractSoilWaterClosure end

closurevar(::SaturationPressureClosure) = auxiliary(:pressure_head, XYZ(), units=u"m", desc="Total hydraulic pressure head in m water displaced at standard pressure")

"""
    $TYPEDSIGNATURES

Computes `pressure_head` ``Ψ = ψm + ψz + ψh`` from the current `saturation_water_ice` state.
"""
function closure!(
    state, grid,
    hydrology::SoilHydrologyRRE{NF, SaturationPressureClosure},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    # apply saturation correction
    launch!(grid, XY, adjust_saturation_profile_kernel!, state, hydrology)
    # update water table
    compute_water_table!(state, grid, hydrology)
    # determine pressure head from saturation
    launch!(grid, XYZ, saturation_to_pressure_kernel!, state, hydrology, strat, bgc)
    return nothing
end

"""
    $TYPEDSIGNATURES

Computes `saturation_water_ice` from the current `pressure_head` state.
"""
function invclosure!(
    state, grid,
    hydrology::SoilHydrologyRRE{NF, SaturationPressureClosure},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    # determine saturation from pressure
    launch!(grid, XYZ, pressure_to_saturation_kernel!, state, hydrology, strat, bgc)
    # apply saturation correction
    launch!(grid, XY, adjust_saturation_profile_kernel!, state, hydrology)
    # update water table
    compute_water_table!(state, grid, hydrology)
    return nothing
end

@kernel function pressure_to_saturation_kernel!(
    state, grid,
    hydrology::SoilHydrologyRRE{NF, SaturationPressureClosure},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    i, j, k = @index(Global, NTuple)
    pressure_to_saturation!(i, j, k, grid, state, hydrology, strat, bgc)
end

@kernel function saturation_to_pressure_kernel!(
    state, grid,
    hydrology::SoilHydrologyRRE{NF, SaturationPressureClosure},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    i, j, k = @index(Global, NTuple)
    saturation_to_pressure!(i, j, k, grid, state, hydrology, strat, bgc)
end

@inline function pressure_to_saturation!(
    i, j, k, grid, state,
    hydrology::SoilHydrologyRRE{NF, SaturationPressureClosure},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    fgrid = get_field_grid(grid)
    ψ = state.pressure_head[i, j, k] # assumed given
    # get the elevation (z-coord) of k'th layer and reference (surface)
    z = znode(i, j, k, state.pressure_head)
    z_ref = znode(i, j, fgrid.Nz, state.pressure_head)
    # elevation pressure head
    ψz = z - z_ref
    # compute hydrostatic pressure head assuming impermeable lower boundary
    # TODO: relax this assumption in the future?
    z₀ = state.water_table[i, j, 1]
    ψh = max(0, z₀ - z)
    # remove hydrostatic and elevation components
    ψm = ψ - ψh - ψz
    swrc = get_swrc(hydrology)
    por = porosity(i, j, k, grid, state, strat, bgc)
    vol_water_ice_content = swrc(ψm; θsat=por)
    state.saturation_water_ice[i, j, k] = vol_water_ice_content / por
    return nothing
end

@inline function saturation_to_pressure!(
    i, j, k, grid, state,
    hydrology::SoilHydrologyRRE{NF, SaturationPressureClosure},
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
) where {NF}
    fgrid = get_field_grid(grid)
    sat = state.saturation_water_ice[i, j, k] # assumed given
    # get the elevation (z-coord) of k'th layer and reference (surface)
    z = znode(i, j, k, state.saturation_water_ice)
    z_ref = znode(i, j, fgrid.Nz, state.saturation_water_ice)
    # get inverse of SWRC
    inv_swrc = inv(get_swrc(hydrology))
    por = porosity(i, j, k, grid, state, strat, bgc)
    # compute matric pressure head
    ψm = inv_swrc(sat*por; θsat=por)
    # compute elevation pressure head
    ψz = z - z_ref
    # compute hydrostatic pressure head assuming impermeable lower boundary
    # TODO: can we generalize this for arbitrary lower boundaries?
    z₀ = state.water_table[i, j, 1]
    ψh = max(0, z₀ - z)
    # compute total pressure head as sum of ψh + ψm + ψz
    # note that ψh and ψz will cancel out in the saturated zone
    state.pressure_head[i, j, k] = ψh + ψm + ψz
    return nothing
end
