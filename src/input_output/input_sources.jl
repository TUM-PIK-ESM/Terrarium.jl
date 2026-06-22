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

varpath(::InputSource{NF, name}) where {NF, name} = varpath(name)

"""
    $TYPEDSIGNATURES

Returns the name of the input variable provided by this source, i.e. the last entry of its [`varpath`](@ref).
"""
varname(source::InputSource) = last(varpath(source))

"""
    $TYPEDSIGNATURES

Determine whether the given `source` provides an input variable for the namespace at the given
`scope`, where `scope` is the path of namespace names from the root namespace, i.e. `()` for the
root namespace itself, `(:ns1,)` for its child namespace `ns1`, and so on.
"""
matches_scope(source::InputSource, scope::VarPath) = Base.front(varpath(source)) == scope

"""
    $TYPEDSIGNATURES

Initializes the input source. The `scope` corresponds to the namespace path (relative to the root namespace)
of the input variables in `fields`. Default implementation does nothing.
"""
initialize!(fields, ::InputSource, clock, scope::VarPath = ()) = nothing

"""
    $TYPEDSIGNATURES

Updates the values of input variables stored in `fields` from the given input `source`. The `scope`
corresponds to the namespace path (relative to the root namespace) of the input variables in `fields`.
Default implementation simply returns `nothing`.
"""
update_inputs!(fields, ::InputSource, ::Clock, scope::VarPath = ()) = nothing

"""
Type alias for an `AbstractField` with any X, Y, Z location or grid.
"""
const AnyField{NF} = AbstractField{LX, LY, LZ, G, NF} where {LX, LY, LZ, G}

"""
Container type for wrapping multiple `InputSource`s.
"""
struct InputSources{NF, Sources <: Tuple{Vararg{InputSource{NF}}}} <: InputSource{NF, nothing}
    sources::Sources
end

InputSources(::Type{NF}) where {NF} = InputSources{NF, Tuple{}}(())
InputSources(input::InputSource, others::InputSource...) = InputSources(tuple(input, others...))
InputSources(input::InputSources) = input

variables(sources::InputSources) = tuplejoin(map(variables, sources.sources)...)

varname(::InputSources) = nothing

function initialize!(fields, input::InputSources, clock::Clock, scope::VarPath = ())
    for source in input.sources
        initialize!(fields, source, clock, scope)
    end
    return nothing
end

function update_inputs!(fields, input::InputSources, clock::Clock, scope::VarPath = ())
    for source in input.sources
        update_inputs!(fields, source, clock, scope)
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

function initialize!(fields, source::FieldInputSource, clock::Clock = Clock(time = 0.0), scope::VarPath = ())
    name = varname(source)
    if matches_scope(source, scope) && hasproperty(fields, name)
        field = getproperty(fields, name)
        set!(field, source.field)
    end
    return nothing
end

"""
    $TYPEDSIGNATURES

Create a `FieldInputSource` with the given grid and input variable `fields`. Use it for static input fields.
The `name` can either be a plain `Symbol` or a namespaced path; see [`varpath`](@ref).
"""
function InputSource(grid::AbstractLandGrid{NF}, field::FS; name, units = NoUnits) where {NF, FS <: AnyField{NF}}
    # ensure fields are on the same architecture as the grid
    field = on_architecture(architecture(grid), field)

    # Check that fields match grid
    @assert field.grid == get_field_grid(grid) "Field must have the same grid as the input grid"

    # infer the VarDims and subsequently the Field location from the data dimensions
    dims = Terrarium.vardims(field)

    path = varpath(name)
    return FieldInputSource{NF, path, typeof(dims), typeof(field), typeof(units)}(dims, units, field)
end

"""
    $TYPEDSIGNATURES

Convenience function to create a `FieldInputSource` from a `RingGrids.Field`.
Converts the RingGrids field to an Oceananigans field and then creates the input source.
"""
function InputSource(grid::ColumnRingGrid{NF}, ring_field::RingGrids.AbstractField; name, units = NoUnits) where {NF}
    oceananigans_field = Field(ring_field, grid)
    dims = Terrarium.vardims(oceananigans_field)
    path = varpath(name)
    return FieldInputSource{NF, path, typeof(dims), typeof(oceananigans_field), typeof(units)}(dims, units, oceananigans_field)
end

variables(source::FieldInputSource) = with_scope(varpath(source), input(varname(source), source.dims; units = source.units))

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
    path = varpath(name)
    return FieldTimeSeriesInputSource{NF, path, typeof(dims), typeof(fts), typeof(units)}(dims, units, fts)
end

variables(source::FieldTimeSeriesInputSource) = with_scope(varpath(source), input(varname(source), source.dims; units = source.units))

# to initialize just update the state once at the start time
function initialize!(fields, source::FieldTimeSeriesInputSource, clock::Clock, scope::VarPath = ())
    return update_inputs!(fields, source, clock, scope)
end

function update_inputs!(fields, source::FieldTimeSeriesInputSource, clock::Clock, scope::VarPath = ())
    name = varname(source)
    if matches_scope(source, scope) && hasproperty(fields, name)
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
