# Boundary conditions interface

abstract type AbstractBoundaryConditions end

"""
    get_field_boundary_conditions(bcs::AbstractBoundaryConditions, var::AbstractVariable)
"""
get_field_boundary_conditions(::AbstractBoundaryConditions, ::AbstractVariable) = nothing

struct PrescribedFluxes <: AbstractBoundaryConditions end

@kwdef struct FieldBoundaryConditions{BCS<:NamedTuple} <: AbstractBoundaryConditions
    var_bcs::BCS = (;)
end


