# Convenience dispatches for Oceananigans.launch!
function launch!(grid::AbstractLandGrid, workspec::Symbol, kernel::Function, args...; kwargs...)
    fgrid = get_field_grid(grid)
    return launch!(fgrid.architecture, fgrid, workspec, kernel, args...; kwargs...)
end

function launch!(state, grid::AbstractLandGrid, workspec::Symbol, kernel::Function, args...; kwargs...)
    fgrid = get_field_grid(grid)
    return launch!(fgrid.architecture, fgrid, workspec, kernel, state, grid, args...; kwargs...)
end

"""
Returns the appropriate workspec for the given `AbstractField` or based on the given
field locations.
"""
workspec(::AbstractField{LX, LY, LZ}) where {LX, LY, LZ} = workspec(LX(), LY(), LZ())
workspec(::Center, ::Center, ::Nothing) = :xy
workspec(::Center, ::Center, ::Center) = :xyz

# Helper functions for checking if a `RingGrids` or `Oceananigans` `Field` matches the given grid
field_matches_grid(field, grid) = field.grid == grid

function assert_field_matches_grid(field::Union{RingGrids.AbstractField, AbstractField}, grid)
    return @assert field_matches_grid(field, grid) "Field grid $(typeof(field.grid)) does not match $(typeof(grid))"
end

const RingGridOrField = Union{RingGrids.AbstractGrid, RingGrids.AbstractField}

on_architecture(::GPU, obj::RingGridOrField) = RingGrids.Architectures.on_architecture(RingGrids.Architectures.GPU(), obj)
on_architecture(::CPU, obj::RingGridOrField) = RingGrids.Architectures.on_architecture(RingGrids.Architectures.CPU(), obj)

# Field construction

"""
    Field(
        grid::AbstractLandGrid,
        dims::VarDims,
        boundary_conditions = nothing,
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
        boundary_conditions = nothing,
        args...;
        kwargs...
    )
    # infer the location of the Field on the Oceananigans grid from `dims`
    loc = location(dims)
    FT = Field{map(typeof, loc)...}
    # Specify BCs if defined
    field = if isa(boundary_conditions, FieldBoundaryConditions)
        FT(get_field_grid(grid), args...; boundary_conditions, kwargs...)
    elseif isa(boundary_conditions, NamedTuple)
        # assume that named tuple corresponds to FieldBoundaryConditions positions
        field_bcs = FieldBoundaryConditions(get_field_grid(grid), (Center(), Center(), nothing); boundary_conditions...)
        FT(get_field_grid(grid), args...; boundary_conditions = field_bcs, kwargs...)
    else
        FT(get_field_grid(grid), args...; kwargs...)
    end
    return field
end

"""
    FieldTimeSeries(
        grid::AbstractLandGrid,
        dims::VarDims,
        times=eltype(grid)[]
    )

Construct a `FieldTimeSeries` on the given land `grid` with the given `dims` and `times`.
"""
function FieldTimeSeries(
        grid::AbstractLandGrid,
        dims::VarDims,
        times = eltype(grid)[]
    )
    loc = location(dims)
    return FieldTimeSeries(loc, get_field_grid(grid), times)
end

# Custom operators (consider moving to separate folder/file once there are more?)

"""
    min_zᵃᵃᶠ(i, j, k, grid, x)
    min_zᵃᵃᶠ(i, j, k, grid, f, args...)

Computes the field or function at the vertical (z-axis) face by taking the `min` of the two adjacent vertical layers.
"""
@inline min_zᵃᵃᶠ(i, j, k, grid, c) = @inbounds min(c[i, j, k], c[i, j, k - 1])
@inline min_zᵃᵃᶠ(i, j, k, grid, f, args...) = @inbounds min(f(i, j, k, grid, args...), f(i, j, k - 1, grid, args...))
