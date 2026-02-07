# Saturation ↔ Hydraulic (pressure) head

"""
    $TYPEDEF

Represents a closure relating saturation of water/ice in soil pores to a corresponding
pressure (or hydraulic) head. Note that here "pressure head" is defined to be synonymous
with hydraulic head, i.e. including all both elevation and hydrostatic pressure contributions.
This relation is typically described by soil property-dependent *soil-water retention curve* (SWRC)
which is here defined in implementations of [`AbstractSoilHydraulics`](@ref).
"""
@kwdef struct SaturationPressureClosure <: AbstractSoilWaterClosure end

variables(::SaturationPressureClosure) = (
    auxiliary(:pressure_head, XYZ(), units=u"m", desc="Total hydraulic pressure head in m water displaced at standard pressure"),
)

@propagate_inbounds function pressure_to_saturation!(
    saturation_water_ice, i, j, k, grid, fields,
    closure::SaturationPressureClosure,
    hydrology::SoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    fgrid = get_field_grid(grid)
    ψ = fields.pressure_head[i, j, k] # assumed given
    # get the elevation (z-coord) of k'th layer and reference (surface)
    z = znode(i, j, k, fgrid, Center(), Center(), Center())
    # TODO: we need a more user friendly interface for this...
    z_ref = znode(i, j, fgrid.Nz+1, fgrid, Center(), Center(), Face())
    # elevation pressure head
    ψz = z - z_ref
    # compute hydrostatic pressure head assuming impermeable lower boundary
    # TODO: relax this assumption in the future?
    z₀ = fields.water_table[i, j, 1]
    ψh = max(0, z₀ - z)
    # remove hydrostatic and elevation components
    ψm = ψ - ψh - ψz
    swrc = get_swrc(hydrology)
    por = porosity(i, j, k, grid, fields, strat, bgc)
    vol_water_ice_content = swrc(ψm; θsat=por)
    saturation_water_ice[i, j, k] = vol_water_ice_content / por
    return nothing
end

@propagate_inbounds function saturation_to_pressure!(
    pressure_head, i, j, k, grid, fields,
    closure::SaturationPressureClosure,
    hydrology::SoilHydrology,
    strat::AbstractStratigraphy,
    bgc::AbstractSoilBiogeochemistry
)
    fgrid = get_field_grid(grid)
    sat = fields.saturation_water_ice[i, j, k] # assumed given
    # get the elevation (z-coord) of k'th layer and reference (surface)
    z = znode(i, j, k, fgrid, Center(), Center(), Center())
    z_ref = znode(i, j, fgrid.Nz+1, fgrid, Center(), Center(), Face())
    # get inverse of SWRC
    inv_swrc = inv(get_swrc(hydrology))
    por = porosity(i, j, k, grid, fields, strat, bgc)
    # compute matric pressure head
    ψm = inv_swrc(sat*por; θsat=por)
    # compute elevation pressure head
    ψz = z - z_ref
    # compute hydrostatic pressure head assuming impermeable lower boundary
    # TODO: can we generalize this for arbitrary lower boundaries?
    z₀ = fields.water_table[i, j, 1]
    ψh = max(0, z₀ - z)
    # compute total pressure head as sum of ψh + ψm + ψz
    # note that ψh and ψz will cancel out in the saturated zone
    pressure_head[i, j, k] = ψh + ψm + ψz
    return nothing
end

# Kernels

@kernel inbounds=true function pressure_to_saturation_kernel!(
    out, grid, fields,
    closure::SaturationPressureClosure,
    args...
)
    i, j, k = @index(Global, NTuple)
    pressure_to_saturation!(out.saturation_water_ice, i, j, k, grid, fields, closure, args...)
end

@kernel inbounds=true function saturation_to_pressure_kernel!(
    out, grid, fields,
    closure::SaturationPressureClosure,
    args...
)
    i, j, k = @index(Global, NTuple)
    saturation_to_pressure!(out.pressure_head, i, j, k, grid, fields, closure, args...)
end
