"""
    $TYPEDEF

Base type for input data sources. Implementations of `InputSource` are free to load data
from any arbitrary backend. They expect an `initialize!(fields, ::InputSource)` that is called once at model initialization and an 
`update_inputs!(fields, ::InputSource, ::Clock)` method that is called at every time step. Both default to doing nothing. Implementations should
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

Initializes the input source. Default implementation does nothing.
"""
initialize!(fields, ::InputSource, clock) = nothing

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
struct InputSources{Sources <: Tuple{Vararg{InputSource}}}
    sources::Sources
end

InputSources(sources::InputSource...) = InputSources(Tuple(sources))

variables(sources::InputSources) = tuplejoin(map(variables, sources.sources)...)

function initialize!(fields, sources::InputSources, clock::Clock)
    for source in sources.sources
        initialize!(fields, source, clock)
    end
    return
end

function update_inputs!(fields, sources::InputSources, clock::Clock)
    for source in sources.sources
        update_inputs!(fields, source, clock)
    end
    return
end

"""
    $TYPEDEF

Input source that defines `input` state variables with the given names which
can then be directly modified by the user.
"""
struct FieldInputSource{NF, VD <: VarDims, FS <: Tuple{Vararg{AnyField{NF}}}, names} <: InputSource{NF}
    "Variable dimensions"
    dims::VD

    "Named tuple of input fields"
    fields::NamedTuple{names, FS}
end

function initialize!(fields, source::FieldInputSource, clock = nothing)
    for name in keys(source.fields)
        if hasproperty(fields, name)
            field = getproperty(fields, name)
            source_field = source.fields[name]
            set!(field, source_field)
        end
    end
    return
end

"""
    $SIGNATURES

Create a `FieldInputSource` with the given grid and input variable `fields`. Use it for static input fields.
"""
function InputSource(grid::AbstractLandGrid{NF}, fields::NamedTuple{names, FS}) where {NF, names, FS <: Tuple{Vararg{AnyField{NF}}}}
    # ensure fields are on the same architecture as the grid
    fields = map(fields) do field
        on_architecture(architecture(grid), field)
    end

    # Check that fields match grid
    field_grid = get_field_grid(grid)
    @assert all(field.grid == field_grid for field in values(fields)) "All fields must have the same grid as the input source"

    # infer the VarDims and subsequently the Field location from the data dimensions
    dims = Terrarium.vardims(first(values(fields)))

    return FieldInputSource{NF, typeof(dims), typeof(values(fields)), names}(dims, fields)
end

"""
    $SIGNATURES

Convenience function to create a `FieldInputSource` from `RingGrids.Field` objects.
Converts the RingGrids fields to Oceananigans fields and then creates the input source.
"""
function InputSource(grid::ColumnRingGrid{NF}, ring_fields::NamedTuple{names, RF}) where {NF, names, RF <: Tuple{Vararg{RingGrids.AbstractField}}}
    # Convert RingGrids fields to Oceananigans fields
    oceananigans_fields = map(ring_fields) do ring_field
        Field(ring_field, grid)
    end
    # Infer dimensions from converted fields
    dims = Terrarium.vardims(first(values(oceananigans_fields)))
    return FieldInputSource{NF, typeof(dims), typeof(values(oceananigans_fields)), names}(dims, oceananigans_fields)
end

variables(source::FieldInputSource{NF, VD, FS, names}) where {NF, VD, FS, names} = map(name -> input(name, source.dims), names)

"""
Type alias for a `FieldTimeSeries` with any X, Y, Z location or grid.
"""
const AnyFieldTimeSeries{NF} = FieldTimeSeries{LX, LY, LZ, TI, K, I, D, G, NF} where {LX, LY, LZ, TI, K, I, D, G}

"""
    $TYPEDEF

Input source that reads input fields from pre-specified Oceananigans `FieldTimeSeries`.
"""
struct FieldTimeSeriesInputSource{NF, VD <: VarDims, names, FTS <: Tuple{Vararg{AnyFieldTimeSeries{NF}}}} <: InputSource{NF}
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

# to initialize just update the state once at the start time
function initialize!(fields, source::FieldTimeSeriesInputSource, clock::Clock)
    return update_inputs!(fields, source, clock)
end

function update_inputs!(fields, source::FieldTimeSeriesInputSource, clock::Clock)
    for name in keys(source.fts)
        if hasproperty(fields, name)
            field_t = getproperty(fields, name)
            fts = source.fts[name]
            set!(field_t, fts[Time(clock.time)])
        end
    end
    return
end

# Internal helper method to check that all Field dimensions match
function checkdims(; named_fields...)
    @assert length(named_fields) >= 1 "at least one input field must be provided"
    # infer dimensions of all provided fields
    field_dims = map(vardims, values(named_fields))
    @assert length(field_dims) == 1 || foldl(==, field_dims) "all fields must have matching dimensions"
    return first(field_dims)
end
