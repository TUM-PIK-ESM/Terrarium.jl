using Rasters: Raster, RasterStack, Ti, hasdim, dims

"""
Alias for `RasterStack(filepath; kwargs...)`. Loads a raster-like dataset from the given file.
"""
load_dataset(filepath::String; kwargs...) = RasterStack(filepath)

abstract type AbstractInterpolator end
struct NoInterpolation <: AbstractInterpolator end

"""
    $SIGNATURES

Retrieves the current value of the given scalar or array.
"""
current_value(x::Number, t) = x
current_value(A::AbstractVecOrMat, t) = A
current_value(f::Field, t) = f

mutable struct InputProvider{NF, Grid<:AbstractLandGrid{NF}, Interp, SourceTypes}
    "Spatial grid on which all model `Field`s are defined"
    grid::Grid

    "Container for all input data sources"
    sources::Dict{String, SourceTypes}

    "Interpolator for spatial input data"
    interpolator::Interp

    "Current time propagated to time-dependent inputs; this field needs to be updated at each time step!"
    time::NF # TODO: This is a bit weird... maybe just having callers pass in the time through `current_value` would be cleaner?
end

function InputProvider(
    grid::AbstractLandGrid, stacks::RasterStack...;
    time=zero(eltype(grid)), # starting time
    interp=NoInterpolator(), # no interpolation for now
)
    stack = merge(stacks...)
    sources = map(keys(stack), values(stack)) do name, dataset
        string(name) => dataset
    end
    return InputProvider(grid, Dict(sources), interp, time)
end

function get_input(provider::InputProvider, str::String, dims::VarDims=XY())
    source = provider.sources[str]
    return input_from_source(provider.grid, source, dims)
end

# Simple initial implementation of input_from_source for rasters that directly match the target grid
function input_from_source(provider::InputProvider{NF, Grid, NoInterpolation}, data::Raster, dims::VarDims) where {NF, Grid}
    # we assume that inputs are given always on grid cell centers;
    # TODO: do we need support for inputs that live on grid corners?
    inputloc(::XY) = (Center(), Center())
    inputloc(::XYZ) = (Center(), Center(), Center())
    # Check that grid and data sizes match
    @assert size(grid) == size(data) "Dimensions of grid $(size(grid)) do not match the data $(size(data))"
    # If the data has a time dimension, handle it as time-varying
    if hasdim(data, Ti)
        loc = inputloc(dims)
        times = convert_times(collect(dims(data, Ti)))
        fts = FieldTimeSeries(loc, get_field_grid(grid), times)
        # directly copy data from src into the FieldTimeSeries
        for i in eachindex(times)
            set!(fts[i], data[Ti(i)])
        end
        return TimeVaryingInputField(fts, provider)
    # otherwise treat it as a static field
    else
        input_field = CenterField()
        set!(input_field, data)
        return StaticInputField(input_field)
    end
end

# Convert time types to second offsets
# TODO: This won't work if the forcing datasets have different starting times...
convert_times(ts::AbstractVector{<:AbstractFloat}) = ts
convert_times(ts::AbstractVector{<:TimeType}) = Dates.value(Second(ts .- ts[1]))

"""
Base type for input data/fields.
"""
abstract type AbstractInput end

"""
Represents an input `Field` that is static, i.e. does not change with time.
"""
struct StaticInputField{F<:AbstractField} <: AbstractInput
    "FieldTimeSeries for the given time-varying field"
    field::F
end

@inline current_value(input::StaticInputField, t) = input.field

@inline set!(f::Field, input::StaticInputField) = set!(f, input.field)

"""
Wraps an Oceananigans `FieldTimeSeries`
"""
struct TimeVaryingInputField{FTS<:FieldTimeSeries, Provider} <: AbstractInput
    "FieldTimeSeries for the given time-varying field"
    field::FTS
    
    "Back-reference to InputProvider to get current time"
    provider::Provider
end

@inline current_value(input::TimeVaryingInputField, t) = input.field[Time(t)]

@inline set!(f::Field, input::TimeVaryingInputField) = set!(f, current_value(input, input.provider.time))

# The nested input provider is mutable so here we adapt it away by just immediately computing the value of the
# Field at the current time
Adapt.adapt_structure(to, input::TimeVaryingInputField) = adapt(to, current_value(input, input.provider.time))
