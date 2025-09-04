"""
    $TYPEDEF

Container type for holding data `Field`s for input variables. `InputFields` stores
the `Field`s for each variable in separate dictionaries for 2D and 3D variables defined on the given
`grid`. This allows for input variable `Field`s to be lazily allocated as-needed.
"""
struct InputFields{NF, Grid<:AbstractLandGrid{NF}, Input2D<:AbstractField, Input3D<:AbstractField}
    "Grid on which the `Field`s are defined"
    grid::Grid

    "Mapping of 2D (XY) input variable names to their corresponding `Field`s"
    inputs2D::OrderedDict{Symbol, Input2D}

    "Mapping of 3D (XYZ) input variable names to their corresponding `Field`s"
    inputs3D::OrderedDict{Symbol, Input3D}
end

function InputFields(grid::AbstractLandGrid)
    # create dummy fields to get type info
    dummy_field2D = Field(grid, XY())
    dummy_field3D = Field(grid, XYZ())
    return InputFields(
        grid,
        OrderedDict{Symbol, typeof(dummy_field2D)}(),
        OrderedDict{Symbol, typeof(dummy_field3D)}()
    )
end

get_input_fields(inputs::InputFields, ::XY) = inputs.inputs2D
get_input_fields(inputs::InputFields, ::XYZ) = inputs.inputs3D

get_input_field(inputs::InputFields, var::AbstractVariable) = get_input_field(inputs, varname(var), vardims(var))
function get_input_field(inputs::InputFields, name::Symbol, dims::VarDims)
    fields = get_input_fields(inputs, dims)
    return get!(fields, name) do
        # Lazily instantiate field if it doesn't yet exist
        Field(inputs.grid, dims)
    end
end

include("input_sources.jl")

include("input_provider.jl")
