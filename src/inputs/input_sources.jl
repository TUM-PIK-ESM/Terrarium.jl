"""
Base type for input data sources. Implementations of `AbstractInputSource` are free to load data
from any arbitrary backend but are required to implement the `update_inputs!` method.
The type argument `NF` refers to the numeric type of the input data.
"""
abstract type AbstractInputSource{NF} end

"""
    $SIGNATURES

Updates the values of input variables stored in `fields` from the given input `source`.
Default implementation simply returns `nothing`.
"""
update_inputs!(::InputFields, ::AbstractInputSource, ::Clock) = nothing

"""
Type alias for an `AbstractField` with any X, Y, Z location or grid.
"""
const AnyField{NF} = AbstractField{LX, LY, LZ, G, NF} where {LX, LY, LZ, G}

struct FieldInputSource{NF, VD<:VarDims, names, Fields<:Tuple{Vararg{AnyField{NF}}}} <: AbstractInputSource{NF}
    "Variable dimensions"
    dims::VD    
    
    "Named tuple of input `Field`s"
    fields::NamedTuple{names, Fields}
end

function FieldInputSource(dims=nothing; named_fields...)
    dims = validate_dims(dims; named_fields...)
    return FieldInputSource(dims, (; named_fields...))
end

function update_inputs!(inputs::InputFields, source::FieldInputSource, ::Clock)
    fields = get_input_fields(inputs, source.dims)
    # directly set input fields to source fields
    # this avoids memory duplication and allows for direct changes to the souce field without copying
    for name in keys(source.fields)
        fields[name] = source.fields[name]
    end
end

"""
Type alias for a `FieldTimeSeries` with any X, Y, Z location or grid.
"""
const AnyFieldTimeSeries{NF} = FieldTimeSeries{LX, LY, LZ, TI, K, I, D, G, NF} where {LX, LY, LZ, TI, K, I, D, G}

struct FieldTimeSeriesInputSource{NF, VD<:VarDims, names, FTS<:Tuple{Vararg{AnyFieldTimeSeries{NF}}}} <: AbstractInputSource{NF}
    "Variable dimensions"
    dims::VD
    
    "Field time series data"
    fts::NamedTuple{names, FTS}
end

function FieldTimeSeriesInputSource(dims=nothing; named_fts...)
    dims = validate_dims(dims; named_fts...)
    return FieldTimeSeriesInputSource(dims, (; named_fts...))
end

function update_inputs!(inputs::InputFields, source::FieldTimeSeriesInputSource, clock::Clock)
    fields = get_input_fields(inputs, source.dims)
    for name in keys(source.fts)
        field_t = get(fields, name) do
            # lazily instantiate input field if not yet defined
            create_field(inputs.grid, source.dims)
        end
        fts = source.fts[name]
        set!(field_t, fts[clock.time])
    end
end

# Internal helper method to check that all Field dimensions match;
# Currently this does not check that the dimension sizes match... but maybe it should?
function validate_dims(dims::Union{Nothing, VarDims}=nothing; named_fields...)
    @assert length(named_fields) >= 1 "at least one input field must be provided"
    # infer dimensions of all provided fields
    field_dims = map(inferdims, values(named_fields))
    dims = isnothing(dims) ? field_dims[1] : dims
    @assert foldl(==, field_dims, init=dims) "all fields must have matching dimensions $dims"
    return dims
end
