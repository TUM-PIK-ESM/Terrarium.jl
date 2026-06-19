"""
    $TYPEDEF

Represents an objective function for nonlinear solvers. The name `target`
refers to the output `Field` which should updated on each iteration.
"""
struct ObjectiveFunction{target, F, DF}
    func::F

    dfunc::DF

    ObjectiveFunction(func, dfunc, target::Symbol) = new{target, typeof(func), typeof(dfunc)}(func, dfunc)
    ObjectiveFunction(func, target::Symbol) = ObjectiveFunction(func, nothing, target)
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

# Evaluate the fixed-point residual F(x) = g(x) - x, where g(x) is the value of
# the target field produced by `step_func!` when the field is initialized to `x`.
@propagate_inbounds function build_residual(
        out, indices, grid, fields, step_func!, target_field, solver,
        func_args...; func_kwargs...
    )
    function residual!(x)
        target_field[indices...] = x
        step_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
        x_new = target_field[indices...]
        return x_new - x
    end
    return residual!
end

# Solvers

include("fixed_point.jl")
include("newton_raphson.jl")
include("root_solvers.jl")
