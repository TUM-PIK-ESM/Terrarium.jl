module TerrariumRastersExt

using Terrarium
using Terrarium: XY, XYZ # variable dims
using Unitful: NoUnits

using Dates
using DocStringExtensions
using Interpolations
using Oceananigans
using Rasters

import Oceananigans.Architectures: on_architecture

"""
    $TYPEDEF

Input source that reads data from one or mroe (possibly time-varying) `Raster`. This type should generally not be
constructed directly but rather via the `InputSource` constructor interface.
"""
struct RasterInputSource{NF, name, VD, TT, IM <: AbstractVector{Int}, RS <: AbstractRaster{NF}, UT} <: InputSource{NF, name}
    "Variable dimensions"
    dims::VD

    "Physical units"
    units::UT

    "Index map for the grid mask"
    idxmap::IM

    "Reference time"
    reftime::TT

    "Data raster"
    raster::RS
end

"""
    InputSource(data::Raster, grid::ColumnRingGrid; name = data.name)

Creates a new `RasterInputSource` from the given `Raster` data and `grid`.
"""
function Terrarium.InputSource(grid::ColumnRingGrid{NF}, raster::AbstractRaster{NF}; name = raster.name, units = NoUnits, reftime = nothing) where {NF}
    # get indices from grid mask
    idxmap = on_architecture(architecture(grid), findall(Array(grid.mask)))
    # infer the VarDims and subsequently the Field location from the data dimensions
    vd = Terrarium.vardims(raster)
    # infer reference time
    reftime = default_reftime(raster, reftime)
    return RasterInputSource{NF, name, typeof(vd), typeof(reftime), typeof(idxmap), typeof(raster), typeof(units)}(vd, units, idxmap, reftime, raster)
end

Terrarium.variables(source::RasterInputSource{NF, name}) where {NF, name} = (Terrarium.input(name, source.dims; units = source.units),)

function Terrarium.initialize!(fields, source::RasterInputSource{NF, name}, clock::Clock) where {NF, name}
    if hasproperty(fields, name)
        field = getproperty(fields, name)
        timedim = dims(source.raster, Ti)
        current_time = timestamp(source.reftime, clock.time)
        initialize_from_raster!(field, source.raster, source.idxmap, timedim, current_time)
    end
    return
end

# for static rasters initialize once and then don't update anymore
function initialize_from_raster!(field, raster, idxmap, timedim::Nothing, current_time)
    field .= view(raster, idxmap)
    return nothing
end

# for time-varying rasters this just updates once at the start time
initialize_from_raster!(field, raster, idxmap, timedim, current_time) = update_from_raster!(field, raster, idxmap, timedim, current_time)

function Terrarium.update_inputs!(fields, source::RasterInputSource{NF, name}, clock::Clock) where {NF, name}
    if hasproperty(fields, name)
        field = getproperty(fields, name)
        timedim = dims(source.raster, Ti)
        current_time = timestamp(source.reftime, clock.time)
        update_from_raster!(field, source.raster, source.idxmap, timedim, current_time)
    end
    return
end

# For static raster we don't need to update
function update_from_raster!(
        field::Field,
        raster::AbstractRaster,
        idxmap::AbstractVector{Int},
        ::Nothing,
        t
    )
    return nothing
end

# Update rule for dynamic (time-varying) raster inputs
function update_from_raster!(
        field::Field,
        raster::AbstractRaster,
        idxmap::AbstractVector{Int},
        timedim::TimeDim,
        t::TimeType
    ) where {TimeType, TimeDim <: Ti{<:AbstractVector{TimeType}}}
    arch = architecture(idxmap)
    indexes = searchsorted(timedim.val, t)
    left, right = last(indexes), first(indexes)
    return @inbounds if left >= 1 && right <= length(timedim)
        # Linear interpolation between points
        x1 = on_architecture(arch, raster[Ti(left)])[idxmap]
        x2 = on_architecture(arch, raster[Ti(right)])[idxmap]
        t1 = timedim[left]
        t2 = timedim[right]
        Δt = Terrarium.convert_dt(t2 - t1)
        ϵ = Terrarium.convert_dt(t - t1)
        x_interp = Δt > 0 ? x1 + ϵ * (x2 - x1) / Δt : x2
        set!(field, x_interp)
    else
        # Note: this implicitly results in flat extrapolation beyond the bounds of the time axis;
        # We may want to make this configurable in the future.
        set!(field, on_architecture(arch, raster[Ti(min(right, length(timedim)))])[idxmap])
    end
end

# conversions of simulation time to reference time scale
timestamp(reftime::DateTime, time::Number) = reftime + Nanosecond(round(time * 1.0e9))
timestamp(reftime::DateTime, time::DateTime) = time
timestamp(reftime::Number, time::Number) = reftime + time
timestamp(reftime::Nothing, time) = time

# Use specified reference time, if provided (i.e. not nothing)
default_reftime(data::AbstractRaster, reftime) = reftime
# Otherwise, take first value from time dimension
default_reftime(data::AbstractRaster, ::Nothing) = default_reftime(dims(data, Ti))
default_reftime(timedim::Ti) = first(timedim)
# If no time dimension is defined, return nothing
default_reftime(::Nothing) = nothing

# Infer VarDims based on the axes defined in the Raster
Terrarium.vardims(A::AbstractDimArray) = Terrarium.vardims(dims(A, X, Y, Z)...)
Terrarium.vardims(::X, ::Y) = XY()
Terrarium.vardims(::X, ::Y, ::Z) = XYZ()

# Adapt rules for Rasters
on_architecture(to, raster::AbstractRaster) = rebuild(
    raster,
    data = on_architecture(to, raster.data),
    dims = map(d -> rebuild(d, val = on_architecture(to, d.val)), dims(raster))
)

end
