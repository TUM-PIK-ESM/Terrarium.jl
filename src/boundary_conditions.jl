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
        # automatically forward call to PrescribedBC
        compute_auxiliary!(state, model, PrescribedBC(name, bcs[name]))
    end
end

"""
Computes any tendency contributions from the given boundary conditions.
"""
compute_tendencies!(state, model, ::AbstractBoundaryConditions) = nothing
compute_tendencies!(state, model, ::BoundaryCondition) = nothing
function compute_tendencies!(state, model, bcs::FieldBCs{names}) where names
    fastiterate(names) do name
        compute_tendencies!(state, model, PrescribedBC(name, bcs[name]))
    end
end

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
struct PrescribedBC{progvar, BC, FT} <: AbstractBoundaryConditions
    "Boundary condition value"
    value::FT

    PrescribedBC(progvar::Symbol, bctype::BCType, value) = new{progvar, typeof(bctype), typeof(value)}(value)
end

variables(bc::PrescribedBC) = variables(bc.value)

# compute_tendencies! for flux-valued boundary conditions
function compute_tendencies!(state, model, ::PrescribedBC{progvar, <:Flux}) where {progvar}
    tend = getproperty(state, tendencyof(progvar))
    prog = getproprerty(state, progvar)
    arch = architecture(get_grid(model))
    clock = state.clock
    compute_z_bcs!(tend, prog, arch, clock, state)
end

function get_field_boundary_conditions(bc::PrescribedBC{progvar, BC}) where {progvar, BC}
    (; progvar => BoundaryCondition(BC(), bc.value))
end

function get_field_boundary_conditions(bc::PrescribedBC{progvar, BC, <:Input{name}}) where {progvar, name, BC}
    (; progvar => BoundaryCondition(BC(), (x, z, input) -> input, field_dependencies=(name,)))
end

# Convenience aliases for PrescribedBC
NoFlux(progvar::Symbol) = PrescribedBC(progvar, Flux(), nothing)
PrescribedFlux(progvar::Symbol, value; kwargs...) = PrescribedBC(progvar, Flux(), value)
PrescribedValue(progvar::Symbol, value; kwargs...) = PrescribedBC(progvar, Value(), value)
PrescribedGradient(progvar::Symbol, value; kwargs...) = PrescribedBC(progvar, Gradient(), value)

"""
    $TYPEDEF

Represents the boundary conditions applied at the top and bottom of a 1D column model
discrietized along the vertical (depthwise) axis.

Properties:
$TYPEDFIELDS
"""
@kwdef struct ColumnBoundaryConditions{
    TopBC,
    BottomBC
} <: AbstractBoundaryConditions
    "Boundary condition(s) applied at the top of the vertical column."
    top::TopBC = DefaultBoundaryConditions()
    
    "Boundary condition(s) applied at the botom of the vertical column."
    bottom::BottomBC = DefaultBoundaryConditions()
end

function get_field_boundary_conditions(
    bcs::ColumnBoundaryConditions,
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

"""
Alias for `ColumnBoundaryConditions`
"""
const ColumnBCs{Top, Bottom} = ColumnBoundaryConditions{Top, Bottom} where {Top, Bottom}

variables(bcs::ColumnBoundaryConditions) = tuplejoin(variables(bcs.top), variables(bcs.bottom))

function compute_auxiliary!(state, model, bcs::ColumnBoundaryConditions)
    compute_auxiliary!(state, model, bcs.top)
    compute_auxiliary!(state, model, bcs.bottom)
end

"""
    FieldBoundaryConditions(grid::AbstractLandGrid, loc::Tuple; at...)

Creates a regularized `PrescribedBCs` type from the given keyword arugments of Oceananigans `BoundaryCondition`s
with keys corresponding to their positions on the domain (i.e. `top`, `bottom`, etc.), as well as the grid and location `loc`.
The location refers the position on the staggered grid at which the boundary conditions are defined. For 1D (vertical) domains,
this is usually `(Center(), Center(), nothing)`.
"""
function FieldBoundaryConditions(grid::AbstractLandGrid, loc::Tuple; at...)
    field_grid = get_field_grid(grid)
    bcs = map(((k, bc)) -> k => isnothing(bc) ? NoFluxBoundaryCondition() : bc, keys(at), values(at))
    # create the PrescribedBCs type
    field_bcs = FieldBoundaryConditions(field_grid, loc; immersed=DefaultBoundaryCondition(), bcs...)
    # return the regularized boundary conditions
    return regularize_field_boundary_conditions(field_bcs, field_grid, loc)
end
