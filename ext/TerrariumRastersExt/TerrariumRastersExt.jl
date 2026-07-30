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
import RingGrids

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

Creates a new `RasterInputSource` from the given `Raster` data and `grid`. The `name` can either
be a plain `Symbol` or a namespaced path; see [`Terrarium.varpath`](@ref).
"""
function Terrarium.InputSource(grid::ColumnRingGrid{NF}, raster::AbstractRaster{NF}; name = raster.name, units = NoUnits, reftime = nothing) where {NF}
    # get indices from grid mask
    idxmap = on_architecture(architecture(grid), findall(Array(grid.mask)))
    # infer the VarDims and subsequently the Field location from the data dimensions
    vd = Terrarium.vardims(raster)
    # infer reference time
    reftime = default_reftime(raster, reftime)
    path = Terrarium.varpath(name)
    return RasterInputSource{NF, path, typeof(vd), typeof(reftime), typeof(idxmap), typeof(raster), typeof(units)}(vd, units, idxmap, reftime, raster)
end

function Terrarium.load_asset(path, name, grid::RingGrids.AbstractGrid, ::Terrarium.NetCDF, ::Type{NF}; indices = (:, :, :), fill_value = NF(NaN)) where {NF}
    raster = replace(Raster(path, name = name), missing => fill_value)[indices...]
    data = reconcile_latitudes(raster, grid)
    field = RingGrids.Field(grid, size(data)[3:end]...)
    field .= reshape(data, :, size(data)[3:end]...)
    return field
end

Terrarium.variables(source::RasterInputSource) = Terrarium.with_scope(
    Base.front(Terrarium.varpath(source)),
    Terrarium.input(Terrarium.varname(source), source.dims; units = source.units)
) |> tuple

function Terrarium.initialize!(inputs, grid, clock, fields, source::RasterInputSource)
    name = Terrarium.varname(source)
    if hasproperty(inputs, name)
        field = getproperty(inputs, name)
        timedim = dims(source.raster, Ti)
        current_time = timestamp(source.reftime, clock.time)
        initialize_from_raster!(field, source.raster, source.idxmap, timedim, current_time)
    end
    return nothing
end

# for static rasters initialize once and then don't update anymore
function initialize_from_raster!(field, raster, idxmap, timedim::Nothing, current_time)
    field .= view(raster, idxmap)
    return nothing
end

# for time-varying rasters this just updates once at the start time
initialize_from_raster!(field, raster, idxmap, timedim, current_time) = update_from_raster!(field, raster, idxmap, timedim, current_time)

function Terrarium.update_inputs!(inputs, grid, clock, fields, source::RasterInputSource)
    name = Terrarium.varname(source)
    if hasproperty(inputs, name)
        field = getproperty(inputs, name)
        timedim = dims(source.raster, Ti)
        current_time = timestamp(source.reftime, clock.time)
        update_from_raster!(field, source.raster, source.idxmap, timedim, current_time)
    end
    return nothing
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

# Reconcile the latitude dimension (assumed to be the second, following the (lon, lat, ...)
# layout used above) with the target grid. Regular lon-lat grids (e.g. ERA5-Land) place points
# on both poles, whereas RingGrids full grids do not, so such data carries two extra latitude
# rings. Drop the first and last latitude (the poles) when present; otherwise require an exact
# match so mismatched data is not silently truncated.
function reconcile_latitudes(data::AbstractArray, grid::RingGrids.AbstractFullGrid)
    nlat_data = size(data, 2)
    nlat_grid = RingGrids.get_nlat(grid)
    if nlat_data == nlat_grid
        return data
    elseif nlat_data == nlat_grid + 2
        trailing = ntuple(_ -> Colon(), ndims(data) - 2)
        return data[:, 2:(nlat_data - 1), trailing...]
    else
        error("Cannot reconcile source latitude count ($nlat_data) with grid latitude count ($nlat_grid)")
    end
end

end
