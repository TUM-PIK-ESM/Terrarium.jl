"""
    $TYPEDEF

Base type for input data sources. Implementations of `InputSource` are free to load data
from any arbitrary backend but are required to implement the
`update_inputs!(fields, ::InputSource, ::Clock)` method. Implementations should
additionally provide a constructor as a dispatch of `InputSource`.

The type argument `NF` corresponds to the numeric type of the input data.
"""
abstract type InputSource{NF} end

# Default kwarg constructor for convenience
InputSource(; kwargs...) = InputSource(kwargs...)

"""
    $SIGNATURES

Returns a tuple of `Symbol`s corresponding to variable names supported by this `InputSource`.
"""
variables(::InputSource) = ()

"""
    $SIGNATURES

Updates the values of input variables stored in `fields` from the given input `source`.
Default implementation simply returns `nothing`.
"""
update_inputs!(fields, ::InputSource, ::Clock) = nothing

"""
Type alias for an `AbstractField` with any X, Y, Z location or grid.
"""
const AnyField{NF} = AbstractField{LX, LY, LZ, G, NF} where {LX, LY, LZ, G}

"""
Container type for wrapping multiple `InputSource`s.
"""
struct InputSources{Sources<:Tuple{Vararg{InputSource}}}
    sources::Sources
end

InputSources(sources::InputSource...) = InputSources(Tuple(sources))

variables(sources::InputSources) = tuplejoin(map(variables, sources.sources)...)

function update_inputs!(fields, sources::InputSources, clock::Clock)
    for source in sources.sources
        update_inputs!(fields, source, clock)
    end
end

"""
    $TYPEDEF

Input source that defines `input` state variables with the given names which
can then be directly modified by the user.
"""
struct FieldInputSource{NF, VD<:VarDims, names} <: InputSource{NF}
    "Variable dimensions"
    dims::VD    

    FieldInputSource(::Type{NF}, dims::VarDims, names::Symbol...) where {NF} = new{NF, typeof(dims), names}(dims)
end

"""
    InputSource(::Type{NF}, names::Symbol...; dims = XY())

Create a `FieldInputSource` with the given numeric type and input variable `names`.
"""
function InputSource(::Type{NF}, names::Symbol...; dims = XY()) where {NF}
    return FieldInputSource(NF, dims, names...)
end

variables(source::FieldInputSource{NF, VD, names}) where {NF, VD, names} = map(name -> input(name, source.dims), names)

# For single fields; users can write directly into the allocated input variable
update_inputs!(fields, ::FieldInputSource, ::Clock) = nothing

"""
Type alias for a `FieldTimeSeries` with any X, Y, Z location or grid.
"""
const AnyFieldTimeSeries{NF} = FieldTimeSeries{LX, LY, LZ, TI, K, I, D, G, NF} where {LX, LY, LZ, TI, K, I, D, G}

"""
    $TYPEDEF

Input source that reads input fields from pre-specified Oceananigans `FieldTimeSeries`.
"""
struct FieldTimeSeriesInputSource{NF, VD<:VarDims, names, FTS<:Tuple{Vararg{AnyFieldTimeSeries{NF}}}} <: InputSource{NF}
    "Variable dimensions"
    dims::VD
    
    "Field time series data"
    fts::NamedTuple{names, FTS}
end

function InputSource(named_fts::Pair{Symbol, <:FieldTimeSeries}...)
    dims = checkdims(; named_fts...)
    return FieldTimeSeriesInputSource(dims, (; named_fts...))
end

variables(source::FieldTimeSeriesInputSource{NF, VD, names}) where {NF, VD, names} = map(name -> input(name, source.dims), names)

function update_inputs!(fields, source::FieldTimeSeriesInputSource, clock::Clock)
    for name in keys(source.fts)
        if hasproperty(fields, name)
            field_t = getproperty(fields, name)
            fts = source.fts[name]
            set!(field_t, fts[Time(clock.time)])
        end
    end
end

# Internal helper method to check that all Field dimensions match
function checkdims(; named_fields...)
    @assert length(named_fields) >= 1 "at least one input field must be provided"
    # infer dimensions of all provided fields
    field_dims = map(vardims, values(named_fields))
    @assert length(field_dims) == 1 || foldl(==, field_dims) "all fields must have matching dimensions"
    return first(field_dims)
end
