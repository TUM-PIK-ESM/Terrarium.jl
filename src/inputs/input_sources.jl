"""
    $TYPEDEF

Base type for input data sources. Implementations of `InputSource` are free to load data
from any arbitrary backend but are required to implement the
`update_inputs!(::InputFields, ::InputSource, ::Clock)` method. Implementations should
additionally provide a constructor as a dispatch of `InputSource`.

The type argument `NF` corresponds to the numeric type of the input data.
"""
abstract type InputSource{NF} end

"""
    InputSource(; inputs...)

Constructs a subtype of `InputSource`s from named inputs given as
keyword arguments. The return type will depend on the type of the input arguments.
"""
InputSource(; inputs...) = InputSource(inputs...)

"""
    $SIGNATURES

Updates the values of input variables stored in `fields` from the given input `source`.
Default implementation simply returns `nothing`.
"""
update_inputs!(::InputFields, ::InputSource, ::Clock) = nothing

"""
Type alias for an `AbstractField` with any X, Y, Z location or grid.
"""
const AnyField{NF} = AbstractField{LX, LY, LZ, G, NF} where {LX, LY, LZ, G}

"""
    $TYPEDEF

Input source that reads directly from pre-specified Oceananigans `Field`s.
"""
struct FieldInputSource{NF, VD<:VarDims, names, Fields<:Tuple{Vararg{AnyField{NF}}}} <: InputSource{NF}
    "Variable dimensions"
    dims::VD    
    
    "Named tuple of input `Field`s"
    fields::NamedTuple{names, Fields}
end

function InputSource(named_fields::Pair{Symbol, <:Field}...)
    dims = checkdims(; named_fields...)
    return FieldInputSource(dims, (; named_fields...))
end

function update_inputs!(inputs::InputFields, source::FieldInputSource, ::Clock)
    for name in keys(source.fields)
        field = get_input_field(inputs, name, source.dims)
        set!(field, source.fields[name])
    end
end

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

function update_inputs!(inputs::InputFields, source::FieldTimeSeriesInputSource, clock::Clock)
    fields = get_input_fields(inputs, source.dims)
    for name in keys(source.fts)
        field_t = get!(fields, name) do
            # lazily instantiate input field if not yet defined
            Field(inputs.grid, source.dims)
        end
        fts = source.fts[name]
        set!(field_t, fts[Time(clock.time)])
    end
end

# Internal helper method to check that all Field dimensions match
function checkdims(; named_fields...)
    @assert length(named_fields) >= 1 "at least one input field must be provided"
    # infer dimensions of all provided fields
    field_dims = map(inferdims, values(named_fields))
    @assert length(field_dims) == 1 || foldl(==, field_dims) "all fields must have matching dimensions"
    return first(field_dims)
end
