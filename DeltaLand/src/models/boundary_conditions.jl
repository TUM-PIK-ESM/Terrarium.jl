# Boundary conditions interface

abstract type AbstractBoundaryConditions end

"""
    get_field_boundary_conditions(bcs::AbstractBoundaryConditions, var::AbstractVariable)

Retrieve the `Field` boundary conditions for the given variable. Defaults to returning `nothing`.
"""
get_field_boundary_conditions(::AbstractBoundaryConditions, ::AbstractVariable) = nothing

struct PrescribedFluxes <: AbstractBoundaryConditions end

@kwdef struct FieldBoundaryConditions{BCS<:NamedTuple} <: AbstractBoundaryConditions
    var_bcs::BCS = (;)
end

get_field_boundary_conditions(bcs::FieldBoundaryConditions, var::AbstractVariable) = get(bcs.var_bcs, varname(var), nothing)
