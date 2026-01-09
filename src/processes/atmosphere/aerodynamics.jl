"""
    $TYPEDEF

Dummy implementation of the aerodynamic resistance that simply returns a constant value.
"""
@kwdef struct ConstantAerodynamicResistance{NF} <: AbstractAerodynamics
    "Constant aerodynamic resistance [s m⁻¹]"
    rₐ::NF = 50.0
end

ConstantAerodynamicResistance(::Type{NF}; kwargs...) where {NF} = ConstantAerodynamicResistance{NF}(; kwargs...)

@inline aerodynamic_resistance(i, j, state, res::ConstantAerodynamicResistance) = res.rₐ
