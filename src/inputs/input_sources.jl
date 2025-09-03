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
const AnyField{NF} = AbstractField{LX, LY, LZ, G, NF} where {LX, LY, LZ, G, NF}

struct FieldInput{NF, VD<:VarDims, FT<:AnyField{NF}} <: AbstractInputSource{NF}
    "Name of the input variable"
    name::Symbol
    
    "Field data"
    field::FT

    "Variable dimensions"
    dims::VD

    FieldInput(name::Symbol, field::Field, dims::VarDims=inferdims(field)) = new{eltype(field), typeof(dims), typeof(field)}(name, field, dims)
end

function update_inputs!(inputs::InputFields, source::FieldInput, ::Clock)
    fields = get_input_fields(inputs, source.dims)
    # directly set input field to source.field
    # this avoids memory duplication and allows for direct changes to the souce field without copying
    fields[source.name] = source.field
end

"""
Type alias for a `FieldTimeSeries` with any X, Y, Z location or grid.
"""
const AnyFieldTimeSeries{NF} = FieldTimeSeries{LX, LY, LZ, TI, K, I, D, G, NF} where {LX, LY, LZ, TI, K, I, D, G, NF}

struct FieldTimeSeriesInput{NF, VD<:VarDims, FTS<:AnyFieldTimeSeries{NF}} <: AbstractInputSource{NF}
    "Name of the input variable"
    name::Symbol
    
    "Field time series data"
    fts::FTS

    "Variable dimensions"
    dims::VD

    FieldTimeSeriesInput(
        name::Symbol,
        fts::FieldTimeSeries,
        dims::VarDims=inferdims(fts)
    ) = new{eltype(field), typeof(dims), typeof(field)}(name, field, dims)
end

function update_inputs!(inputs::InputFields, source::FieldTimeSeriesInput, clock::Clock)
    fields = get_input_fields(inputs, source.dims)
    field_t = get(fields, source.name) do
        # lazily instantiate input field if not yet defined
        create_field(inputs.grid, source.dims)
    end
    set!(field_t, source.fts[clock.time])
end
