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
compute_auxiliary!(state, model, ::FieldBCs) = nothing

"""
Computes any tendency contributions from the given boundary conditions.
"""
compute_tendencies!(state, model, ::AbstractBoundaryConditions) = nothing
compute_tendencies!(state, model, ::BoundaryCondition) = nothing
compute_tendencies!(state, model, ::FieldBCs) = nothing

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
struct PrescribedBC{progvar, BC<:BoundaryCondition} <: AbstractBoundaryConditions
    "Boundary condition"
    condition::BC

    PrescribedBC(progvar::Symbol, bc::BoundaryCondition) = new{progvar, typeof(bc)}(bc)
end

variables(bc::PrescribedBC) = variables(bc.condition)

compute_auxiliary!(state, model, ::PrescribedBC) = nothing
compute_tendencies!(state, model, ::PrescribedBC) = nothing

# compute_tendencies! for flux-valued boundary conditions
function compute_tendencies!(state, model, ::PrescribedBC{progvar, <:Flux}) where {progvar}
    tend = getproperty(state.tendencies, progvar)
    prog = getproprerty(state, progvar)
    arch = architecture(get_grid(model))
    clock = state.clock
    compute_z_bcs!(tend, prog, arch, clock, state)
end

function get_field_boundary_conditions(bc::PrescribedBC{progvar}, ::AbstractLandGrid) where {progvar}
    (; progvar => bc.condition)
end

# Convenience aliases for PrescribedBC
NoFlux(progvar::Symbol) = PrescribedBC(progvar, NoFluxBoundaryCondition())
PrescribedFlux(progvar::Symbol, value; kwargs...) = PrescribedBC(progvar, FluxBoundaryCondition(value; kwargs...))
PrescribedValue(progvar::Symbol, value; kwargs...) = PrescribedBC(progvar, ValueBoundaryCondition(value; kwargs...))
PrescribedGradient(progvar::Symbol, value; kwargs...) = PrescribedBC(progvar, GradientBoundaryCondition(value; kwargs...))

"""
Implementation of `Oceananigans.BoundaryConditions.getbc` for `Input{name}` placeholders that retrieves the input `Field` from
`state.inputs` and returns the value at the given index.
"""
@inline function getbc(::Input{name, units, <:XY}, i::Integer, j::Integer, grid::OceananigansGrids.AbstractGrid, clock, state) where {name, units}
    input_field = getproperty(state.inputs, name)
    return @inbounds input_field[i, j]
end

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
        bc = map(x -> isnothing(x) ? NoFluxBoundaryCondition() : x, bc)
        FieldBoundaryConditions(get_field_grid(grid), (Center(), Center(), nothing); bc...)
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
