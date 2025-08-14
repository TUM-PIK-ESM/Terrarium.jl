abstract type AbstractBoundaryConditions end

"""
    $SIGNATURES

Returns an `Oceananigans` `FieldBoundaryConditions` type describing the boundary conditions at
each boundary of the `Field` for the given state variable. Defaults to returning `nothing` which
will result in default boundary conditions being assigned.
"""
get_field_boundary_conditions(::AbstractBoundaryConditions, ::AbstractLandGrid, ::AbstractVariable) = nothing

"""
    get_field_boundary_conditions(bcs::NamedTuple, var::AbstractVariable, grid::AbstractLandGrid, loc::Tuple)

Creates a regularized `FieldBoundaryConditions` type from the given named tuple of Oceananigans boundary condition types
with keys corresponding to their positions on the domain (i.e. `top`, `bottom`, etc.), as well as the grid and location `loc`.
The location refers the position on the staggered grid at which the boundary conditions are defined. For 1D (vertical) domains,
this is usually `(Center(), Center(), nothing)`.
"""
function get_field_boundary_conditions(bcs::NamedTuple, grid::AbstractLandGrid, loc::Tuple=(Center(), Center(), nothing))
    field_grid = get_field_grid(grid)
    field_bcs = FieldBoundaryConditions(field_grid, loc; bcs...)
    return regularize_field_boundary_conditions(field_bcs, field_grid, loc)
end

"""
Marker type for using default boundary conditions inferred from the `Field` and grid types.
"""
struct DefaultBoundaryConditions <: AbstractBoundaryConditions end

"""
    $TYPEDEF

Represents a choice of boundary conditions for a specific state variable with the given `name`. Currently only 1D (top and bottom)
bounday conditions are supported.

Properties:
$TYPEDFIELDS
"""
struct VarBoundaryConditions{TopBC, BottomBC} <: AbstractBoundaryConditions
    "Name of the state variable on which these boundary conditions are defined"
    name::Symbol

    "Boundary conditions at the top of the spatial domain"
    top::TopBC

    "Boundary conditions at the bottom of the spatial domain"
    bottom::BottomBC
end

VarBoundaryConditions(name::Symbol; top=NoFluxBoundaryCondition(), bottom=NoFluxBoundaryCondition()) = VarBoundaryConditions(name, top, bottom)

function get_field_boundary_conditions(bc::VarBoundaryConditions, grid::AbstractLandGrid, var::AbstractVariable)
    if varname(var) == bc.name
        return get_field_boundary_conditions((top=bc.top, bottom=bc.bottom), grid)
    else
        return nothing
    end
end

"""
    $TYPEDEF

Container type that bundle multiple `AbstractBoundaryConditions` structs into a single object that can be supplied to a model.
"""
struct BoundaryConditions{names, BCS<:Tuple{Vararg{AbstractBoundaryConditions}}} <: AbstractBoundaryConditions
    bcs::NamedTuple{names, BCS}
end

BoundaryConditions(; bcs...) = BoundaryConditions((; bcs...))
BoundaryConditions(bcs::VarBoundaryConditions...) = BoundaryConditions(; map(bc -> varname(bc) => bc, bcs)...)

Base.propertynames(bcs::BoundaryConditions) = propertynames(getfield(bcs, :bcs))
Base.getproperty(bcs::BoundaryConditions, name::Symbol) = getproperty(getfield(bcs, :bcs), name)

function get_field_boundary_conditions(bcs::BoundaryConditions, grid::AbstractLandGrid, var::AbstractVariable)
    field_bcs = map(getfield(bcs, :bcs)) do bc
        get_field_boundary_conditions(bc, grid, var)
    end
    # get all non-nothing values
    matched_idx = findall(!isnothing, field_bcs)
    if length(matched_idx) > 1
        @warn "Found more than one matching boundary condition type for $(varname(var)); selecting $(typeof(field_bcs[matched_idx[1]]))"
    end
    return length(matched_idx) >= 1 ? field_bcs[matched_idx[1]] : nothing
end
