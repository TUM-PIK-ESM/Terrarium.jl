"""
    $TYPEDEF

Dummy implementation of aerodynamics that simply returns constant values for all drag coefficients.
"""
@kwdef struct ConstantAerodynamics{NF} <: AbstractAerodynamics
    "Drag coefficient for heat transfer"
    Cₕ::NF = 1.2e-3
end

ConstantAerodynamics(::Type{NF}; kwargs...) where {NF} = ConstantAerodynamics{NF}(; kwargs...)

@inline drag_coefficient(i, j, state, grid, aero::ConstantAerodynamics) = aero.Cₕ
