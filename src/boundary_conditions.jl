"""
Alias for union type of `FieldBoundaryConditions` or a named tuple of `BoundaryCondition`s with
keys corresponding to boundary locations (i.e. top, bottom, etc.)
"""
const FieldBC = Union{FieldBoundaryConditions, NamedTuple{locs, <:Tuple{Vararg{BoundaryCondition}}}} where {locs}

"""
Alias for a `NamedTuple` of `FieldBC` types where the keys correspond to field/variable names.
"""
const FieldBCs{names, BCs} = NamedTuple{names, BCs} where {names, BCs <: Tuple{Vararg{FieldBC}}}

"""
    boundary_conditions(bcs::FieldBCs...)

Recursively merge an arbitrary number of field/variable boundary conditions.
"""
boundary_conditions(bcs::FieldBCs...) = merge_recursive(bcs...)

"""
Implementation of `Oceananigans.BoundaryConditions.getbc` for variable placeholders that retrieves the input `Field` from
`state` and returns the value at the given index.
"""
@inline function getbc(::Variable{name}, i::Integer, j::Integer, grid::OceananigansGrids.AbstractGrid, clock, state::StateVariables) where {name}
    field = getproperty(state, name)
    return @inbounds field[i, j]
end
