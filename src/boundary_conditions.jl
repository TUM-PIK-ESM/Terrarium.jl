abstract type AbstractBoundaryConditions end

"""
Alias for a `NamedTuple` of `BoundaryCondition` types.
"""
const FieldBCs{names, BCs} = NamedTuple{names, BCs} where {names, BCs<:Tuple{Vararg{BoundaryCondition}}}

"""
    $SIGNATURES

Constructs a suitable Oceananigans `BoundaryCondition` for the given state variable `var` on `grid`.
If `bcs` is a `NamedTuple`, it is assumed that the keys correspond to variable names and this method
is invoked recursively on the entry matching the name of `var`, if it exists. Otherwise, `nothing` is
returned.
"""
get_field_boundary_conditions(::AbstractBoundaryConditions, grid::AbstractLandGrid) = (;)
get_field_boundary_conditions(bc::BoundaryCondition, ::AbstractLandGrid) = bc
get_field_boundary_conditions(bcs::FieldBCs, ::AbstractLandGrid) = bcs

"""
Like models/processes, boundary conditions can define state variables which may be computed from
other state variables or from input data in `compute_auxiliary!`.
"""
variables(::AbstractBoundaryConditions) = ()
variables(bc::BoundaryCondition) = variables(bc.condition)
variables(bc::Union{ContinuousBoundaryFunction, DiscreteBoundaryFunction}) = variables(bc.func)
variables(bcs::FieldBCs) = tuplejoin(map(variables, bcs)...)

"""
Updates any state variables associated with the given boundary conditions.
"""
compute_auxiliary!(state, model, ::AbstractBoundaryConditions) = nothing
compute_auxiliary!(state, model, ::BoundaryCondition) = nothing
function compute_auxiliary!(state, model, bcs::FieldBCs{names}) where names
    fastiterate(names) do name
        # automatically forward call to FieldBoundaryCondition
        compute_auxiliary!(state, model, FieldBoundaryCondition(name, bcs[name]))
    end
end

"""
Computes any tendency contributions from the given boundary conditions.
"""
compute_tendencies!(state, model, ::AbstractBoundaryConditions) = nothing
compute_tendencies!(state, model, ::BoundaryCondition) = nothing
function compute_tendencies!(state, model, bcs::FieldBCs{names}) where names
    fastiterate(names) do name
        compute_tendencies!(state, model, FieldBoundaryCondition(name, bcs[name]))
    end
end

"""
    $SIGNATURES

Computes the boundary tendency for the grid cell at `loc`; zero for all other grid cells.
"""
boundary_tendency(i, j, k, grid, loc, state, bc::AbstractBoundaryConditions, args...) = zero(eltype(grid))

"""
Represents default boundary conditions for all state variables. This typically
implies zero-flux boundary conditions on bounded domains and periodic boundary
conditions on periodic domains.
"""
struct DefaultBoundaryConditions <: AbstractBoundaryConditions end

"""
    $TYPEDEF

Container type for an Oceananigans `BoundaryCondition` applied to prognostic (or closure) variable `progvar`.

Properties:
$TYPEDFIELDS
"""
struct FieldBoundaryCondition{progvar, BC} <: AbstractBoundaryConditions
    "Field boundary condition type"
    bc::BC

    FieldBoundaryCondition(progvar::Symbol, bc::BoundaryCondition) = new{progvar, typeof(bc)}(bc)
end

variables(bc::FieldBoundaryCondition) = variables(bc.bc)

function compute_auxiliary!(state, model, bc::FieldBoundaryCondition)
    compute_auxiliary!(state, model, bc.bc)
end

function compute_tendencies!(state, model, bc::FieldBoundaryCondition)
    compute_tendencies!(state, model, bc.bc)
end

# compute_tendencies! for flux-valued boundary conditions
function compute_tendencies!(state, model, ::FieldBoundaryCondition{progvar, <:BoundaryCondition{Flux}}) where {progvar}
    tend = getproperty(state, tendencyof(progvar))
    prog = getproprerty(state, progvar)
    arch = architecture(get_grid(model))
    clock = state.clock
    compute_z_bcs!(tend, prog, arch, clock, state)
end

get_field_boundary_conditions(bc::FieldBoundaryCondition{progvar}) where {progvar} = (; progvar => bc.bc)

# Convenience aliases for FieldBoundaryCondition
PrescribedFlux(progvar::Symbol, value; kwargs...) = FieldBoundaryCondition(progvar, FluxBoundaryCondition(value; kwargs...))
PrescribedValue(progvar::Symbol, value; kwargs...) = FieldBoundaryCondition(progvar, ValueBoundaryCondition(value; kwargs...))
PrescribedGradient(progvar::Symbol, value); kwargs... = FieldBoundaryCondition(progvar, GradientBoundaryCondition(value; kwargs...))

"""
    $TYPEDEF

Represents the boundary conditions applied at the top and bottom of a 1D column model
discrietized along the vertical (depthwise) axis.

Properties:
$TYPEDFIELDS
"""
@kwdef struct VerticalBoundaryConditions{
    TopBC,
    BottomBC
} <: AbstractBoundaryConditions
    "Boundary condition(s) applied at the top of the vertical column."
    top::TopBC = DefaultBoundaryConditions()
    
    "Boundary condition(s) applied at the botom of the vertical column."
    bottom::BottomBC = DefaultBoundaryConditions()
end

function get_field_boundary_conditions(
    bcs::VerticalBoundaryConditions,
    grid::AbstractLandGrid,
)
    top = get_field_boundary_conditions(bcs.top, grid)
    bottom = get_field_boundary_conditions(bcs.bottom, grid)
    # invert the nested structure of the top/bottom named tuples;
    # this way we have a single named tuple where each entry has `top` and `bottom` as properties
    var_bcs = merge_recursive(map(bc -> (top=bc,), top), map(bc -> (bottom=bc,), bottom))
    return map(var_bcs) do bc
        FieldBoundaryConditions(grid, (Center(), Center(), nothing); bc...)
    end
end

variables(bcs::VerticalBoundaryConditions) = tuplejoin(variables(bcs.top), variables(bcs.bottom))

function compute_auxiliary!(state, model, bcs::VerticalBoundaryConditions)
    compute_auxiliary!(state, model, bcs.top)
    compute_auxiliary!(state, model, bcs.bottom)
end

"""
    FieldBoundaryConditions(grid::AbstractLandGrid, loc::Tuple; at...)

Creates a regularized `FieldBoundaryConditions` type from the given keyword arugments of Oceananigans `BoundaryCondition`s
with keys corresponding to their positions on the domain (i.e. `top`, `bottom`, etc.), as well as the grid and location `loc`.
The location refers the position on the staggered grid at which the boundary conditions are defined. For 1D (vertical) domains,
this is usually `(Center(), Center(), nothing)`.
"""
function FieldBoundaryConditions(grid::AbstractLandGrid, loc::Tuple; at...)
    field_grid = get_field_grid(grid)
    bcs = map(((k, bc)) -> k => isnothing(bc) ? NoFluxBoundaryCondition() : bc, keys(at), values(at))
    # create the FieldBoundaryConditions type
    field_bcs = FieldBoundaryConditions(field_grid, loc; immersed=DefaultBoundaryCondition(), bcs...)
    # return the regularized boundary conditions
    return regularize_field_boundary_conditions(field_bcs, field_grid, loc)
end
