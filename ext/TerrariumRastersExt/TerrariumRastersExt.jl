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

Input source that reads data from a time-varying `Raster`, linearly interpolating in time between records.
When `cycle = true`, the time axis is treated as periodic (e.g. an annual climatology) and repeated over the
whole simulation. This type should generally not be constructed directly but rather via the `InputSource`
constructor interface; static (time-invariant) rasters instead yield a plain `FieldInputSource`.
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

    "Whether to treat the time axis as periodic (cyclic) and repeat it over the simulation"
    cycle::Bool
end

"""
    InputSource(grid::ColumnRingGrid, raster::AbstractRaster; name = raster.name, units, reftime, cycle)

Create an `InputSource` for the given `Raster` data and `grid`. A time-invariant raster (no `Ti` dimension)
is regridded once and returned as a plain `Terrarium.FieldInputSource`; a time-varying raster is returned as
a `RasterInputSource` that linearly interpolates in time, optionally cycling the time axis (`cycle = true`).
The `name` can either be a plain `Symbol` or a namespaced path; see [`Terrarium.varpath`](@ref).
"""
function Terrarium.InputSource(grid::ColumnRingGrid{NF}, raster::AbstractRaster{NF}; name = raster.name, units = NoUnits, reftime = nothing, cycle = false) where {NF}
    return raster_input_source(grid, raster, dims(raster, Ti); name, units, reftime, cycle)
end

# Static (time-invariant) rasters reduce to a plain `FieldInputSource`: regrid once into a `RingGrids.Field`
# and defer to the core input source, eliminating the redundant custom source for static data.
function raster_input_source(grid::ColumnRingGrid{NF}, raster::AbstractRaster{NF}, ::Nothing; name, units, reftime, cycle) where {NF}
    idxmap = findall(Array(grid.mask))
    ring_field = similar(grid.mask, NF)
    fill!(parent(ring_field), zero(NF))
    # scatter the masked (active) cells; `Field(ring_field, grid)` gathers exactly these back out
    @inbounds parent(ring_field)[idxmap] .= view(raster, idxmap)
    return Terrarium.InputSource(grid, ring_field; name, units)
end

# Time-varying rasters retain the custom source with lazy time interpolation (and optional cycling).
function raster_input_source(grid::ColumnRingGrid{NF}, raster::AbstractRaster{NF}, ::Rasters.Dimension; name, units, reftime, cycle) where {NF}
    # get indices from grid mask
    idxmap = on_architecture(architecture(grid), findall(Array(grid.mask)))
    # infer the VarDims and subsequently the Field location from the data dimensions
    vd = Terrarium.vardims(raster)
    # infer reference time
    reftime = default_reftime(raster, reftime)
    path = Terrarium.varpath(name)
    return RasterInputSource{NF, path, typeof(vd), typeof(reftime), typeof(idxmap), typeof(raster), typeof(units)}(vd, units, idxmap, reftime, raster, cycle)
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

# A time-varying source is initialized by writing its value at the start time.
Terrarium.initialize!(inputs, grid, clock, fields, source::RasterInputSource) =
    Terrarium.update_inputs!(inputs, grid, clock, fields, source)

function Terrarium.update_inputs!(inputs, grid, clock, fields, source::RasterInputSource)
    name = Terrarium.varname(source)
    if hasproperty(inputs, name)
        field = getproperty(inputs, name)
        timedim = dims(source.raster, Ti)
        current_time = timestamp(source.reftime, clock.time)
        update_from_raster!(field, source.raster, source.idxmap, timedim, current_time, source.cycle)
    end
    return nothing
end

# Dispatch on the `cycle` flag: cyclic sources wrap the time axis, otherwise use the standard interpolation.
function update_from_raster!(field, raster, idxmap, timedim, t, cycle::Bool)
    return cycle ? update_from_raster_cyclic!(field, raster, idxmap, timedim, t) :
        update_from_raster!(field, raster, idxmap, timedim, t)
end

# Cyclic (periodic) time axis: map the requested time into one period and interpolate, wrapping the
# last record back to the first across the period boundary so the series repeats without a phase drift.
function update_from_raster_cyclic!(field, raster, idxmap, timedim::TimeDim, t) where {TimeDim <: Ti}
    tvals = timedim.val
    N = length(tvals)
    N == 1 && return update_from_raster!(field, raster, idxmap, timedim, t)
    t0 = first(tvals)
    span = Terrarium.convert_dt(last(tvals) - t0)
    Δ = span / (N - 1)          # mean record spacing; closes the loop so the period is a full cycle
    period = span + Δ
    ϕ = mod(Terrarium.convert_dt(t - t0), period)
    if ϕ > span
        # within the wrap gap between the last and first records
        arch = architecture(idxmap)
        w = (ϕ - span) / Δ
        xN = on_architecture(arch, raster[Ti(N)])[idxmap]
        x1 = on_architecture(arch, raster[Ti(1)])[idxmap]
        return set!(field, xN .+ w .* (x1 .- xN))
    else
        # within the sampled range: delegate to the standard interpolation at the wrapped time
        return update_from_raster!(field, raster, idxmap, timedim, add_elapsed(t0, ϕ))
    end
end

# Add an elapsed number of seconds `ϕ` to a reference time `t0` (DateTime or plain number).
add_elapsed(t0::TimeType, ϕ) where {TimeType} = t0 + Nanosecond(round(Int, ϕ * 1.0e9))
add_elapsed(t0::Number, ϕ) = t0 + ϕ

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
