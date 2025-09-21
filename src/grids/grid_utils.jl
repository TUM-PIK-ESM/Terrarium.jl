# Convenience dispatch for Oceananigans.launch!
function launch!(grid::AbstractLandGrid, workdims::Symbol, args...; kwargs...)
    _grid = get_field_grid(grid)
    launch!(_grid.architecture, _grid, workdims, args...; kwargs...)
end

# Field construction

"""
    Field(
        grid::AbstractLandGrid,
        dims::VarDims,
        boundary_conditions=nothing,
        args...;
        kwargs...
    )

Auxiliary constructor for an Oceananigans `Field` on `grid` with the given Terrarium variable `dims` and boundary conditions.
Additional arguments are passed direclty to the `Field` constructor. The location of the `Field`
is determined by `VarDims` defined on `var`.
"""
function Field(
    grid::AbstractLandGrid,
    dims::VarDims,
    boundary_conditions=nothing,
    args...;
    kwargs...
)
    # infer the location of the Field on the FVM grid and specify its type
    loc = inferloc(dims)
    FT = Field{map(typeof, loc)...}
    # Specify BCs if defined
    field = if !isnothing(boundary_conditions)
        FT(get_field_grid(grid), args...; boundary_conditions, kwargs...)
    else
        FT(get_field_grid(grid), args...; kwargs...)
    end
    return field
end

function FieldTimeSeries(
    grid::AbstractLandGrid,
    dims::VarDims,
    times=eltype(grid)[]
)
    loc = inferloc(dims)
    return FieldTimeSeries(loc, get_field_grid(grid), times)
end
