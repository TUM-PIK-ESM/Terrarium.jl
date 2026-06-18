"""
    $TYPEDEF

Represents an objective function for nonlinear solvers. The name `target`
refers to the output `Field` which should updated on each iteration.
"""
struct ObjectiveFunction{target, F}
    func::F

    ObjectiveFunction(func, target::Symbol) = new{target, typeof(func)}(func)
end

@inline (func::ObjectiveFunction)(args...) = func.func(args...)

"""
    $TYPEDEF

Simple relaxation scheme for iterative solvers.
"""
@kwdef struct RelaxationFactor{NF}
    "Relaxation factor"
    factor::NF = 0.5

    "Upper limit for relaxed updates"
    upper_limit::NF = oftype(factor, Inf)

    "Lower limit for relaxed updates"
    lower_limit::NF = oftype(factor, -Inf)
end

"""
    $TYPEDSIGNATURES

Apply the under-relaxed fixed-point update to the skin temperature. 
"""
@inline function relaxed_update(relax::RelaxationFactor{NF}, x_target::NF, x_old::NF) where {NF}
    ω = relax.factor
    x_new = (1 - ω) * x_old + ω * x_target
    return clamp(x_new, relax.lower_limit, relax.upper_limit)
end

relaxed_update(::Nothing, x_target, x_old) = x_target

# Interface

function solve! end

# Solvers

include("fixed_point.jl")
include("newton_raphson.jl")
