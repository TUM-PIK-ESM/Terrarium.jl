"""
    $TYPEDEF

Base type for input data sources. Implementations of `InputSource` are free to load data
from any arbitrary backend. They expect an `initialize!(fields, ::InputSource)` that is called once at model initialization and an 
`update_inputs!(fields, ::InputSource, ::Clock)` method that is called at every time step. Both default to doing nothing. Implementations should
additionally provide a constructor as a dispatch of `InputSource`.

The type argument `NF` corresponds to the numeric type of the input data, `name` to its name that's also used in its `variables` definition. 
"""
abstract type InputSource{NF, name} end

# Default kwarg constructor for convenience
InputSource(; kwargs...) = InputSource(kwargs...)

"""
    $TYPEDSIGNATURES

Returns a tuple of `Symbol`s corresponding to variable names supported by this `InputSource`.
"""
variables(::InputSource) = ()

"""
    $TYPEDSIGNATURES

Returns the name of the input source.
"""
varname(::InputSource{NF, name}) where {NF, name} = name

"""
    $TYPEDSIGNATURES

Initializes the input source. Default implementation does nothing.
"""
initialize!(fields, ::InputSource, clock) = nothing

"""
    $TYPEDSIGNATURES

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
    return nothing
end

function update_inputs!(fields, sources::InputSources, clock::Clock)
    for source in sources.sources
        update_inputs!(fields, source, clock)
    end
    return nothing
end

"""
    $TYPEDEF

Input source that defines `input` state variables with the given names which
can then be directly modified by the user.
"""
struct FieldInputSource{NF, name, VD <: VarDims, FS <: AnyField{NF}, UT} <: InputSource{NF, name}
    "Variable dimensions"
    dims::VD

    "Physical units"
    units::UT

    "Input field"
    field::FS
end

function initialize!(fields, source::FieldInputSource{NF, name}, clock = nothing) where {NF, name}
    if hasproperty(fields, name)
        field = getproperty(fields, name)
        set!(field, source.field)
    end
    return nothing
end

"""
    $TYPEDSIGNATURES

Create a `FieldInputSource` with the given grid and input variable `fields`. Use it for static input fields.
"""
function InputSource(grid::AbstractLandGrid{NF}, field::FS; name, units = NoUnits) where {NF, FS <: AnyField{NF}}
    # ensure fields are on the same architecture as the grid
    field = on_architecture(architecture(grid), field)

    # Check that fields match grid
    @assert field.grid == get_field_grid(grid) "Field must have the same grid as the input grid"

    # infer the VarDims and subsequently the Field location from the data dimensions
    dims = Terrarium.vardims(field)

    return FieldInputSource{NF, name, typeof(dims), typeof(field), typeof(units)}(dims, units, field)
end

"""
    $TYPEDSIGNATURES

Convenience function to create a `FieldInputSource` from a `RingGrids.Field`.
Converts the RingGrids field to an Oceananigans field and then creates the input source.
"""
function InputSource(grid::ColumnRingGrid{NF}, ring_field::RingGrids.AbstractField; name, units = NoUnits) where {NF}
    oceananigans_field = Field(ring_field, grid)
    dims = Terrarium.vardims(oceananigans_field)
    return FieldInputSource{NF, name, typeof(dims), typeof(oceananigans_field), typeof(units)}(dims, units, oceananigans_field)
end

variables(source::FieldInputSource{NF, name}) where {NF, name} = (input(name, source.dims; units = source.units),)

"""
Type alias for a `FieldTimeSeries` with any X, Y, Z location or grid.
"""
const AnyFieldTimeSeries{NF} = FieldTimeSeries{LX, LY, LZ, TI, K, I, D, G, NF} where {LX, LY, LZ, TI, K, I, D, G}

"""
    $TYPEDEF

Input source that reads input fields from pre-specified Oceananigans `FieldTimeSeries`.
"""
struct FieldTimeSeriesInputSource{NF, name, VD <: VarDims, FTS <: AnyFieldTimeSeries{NF}, UT} <: InputSource{NF, name}
    "Variable dimensions"
    dims::VD

    "Physical units"
    units::UT

    "Field time series data"
    fts::FTS
end

function InputSource(fts::AnyFieldTimeSeries{NF}; name, units = NoUnits) where {NF}
    dims = vardims(fts)
    return FieldTimeSeriesInputSource{NF, name, typeof(dims), typeof(fts), typeof(units)}(dims, units, fts)
end

variables(source::FieldTimeSeriesInputSource{NF, name}) where {NF, name} = (input(name, source.dims; units = source.units),)

# to initialize just update the state once at the start time
function initialize!(fields, source::FieldTimeSeriesInputSource, clock::Clock)
    return update_inputs!(fields, source, clock)
end

function update_inputs!(fields, source::FieldTimeSeriesInputSource{NF, name}, clock::Clock) where {NF, name}
    if hasproperty(fields, name)
        field_t = getproperty(fields, name)
        set!(field_t, source.fts[Time(clock.time)])
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
