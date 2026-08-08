"""
    $TYPEDEF

Fixed-point (Picard) iteration solver. Repeatedly applies the update `g(x) = x - F(x)`, where
`F` is the residual returned by the [`ObjectiveFunction`](@ref), optionally under-relaxed via a
[`RelaxationFactor`](@ref), until the change in the iterate falls below `tolerance` or
`max_iterations` is reached.

Properties:
$FIELDS
"""
struct FixedPointSolver{NF, R}
    "Numerical tolerance of the fixed point iteration"
    tolerance::NF

    "Relaxation scheme"
    relax::R

    "Maximum number of iterations to run"
    max_iterations::Int
end

function FixedPointSolver(
        ::Type{NF};
        tolerance::NF = sqrt(eps(NF)),
        relax::R = RelaxationFactor(factor = NF(0.5)),
        max_iterations::Int = 100
    ) where {NF, R}
    return FixedPointSolver{NF, R}(tolerance, relax, max_iterations)
end

@propagate_inbounds function solve!(
        out, indices::NTuple{N, Integer}, grid, fields,
        objective_func!::ObjectiveFunction{target},
        solver::FixedPointSolver{NF},
        func_args...;
        func_kwargs...
    ) where {NF, N, target}
    target_field = getproperty(out, target)
    x₀ = target_field[field_indices...]
    # The objective returns the fixed-point residual F(x) = x - g(x), so the
    # updated iterate is g(x) = x - F(x). Take the first (un-relaxed) step.
    residual = objective_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
    x₁ = x₀ - residual
    target_field[field_indices(indices)...] = x₁
    iteration = 1
    while abs(x₁ - x₀) > solver.tolerance && iteration <= solver.max_iterations
        x₀ = x₁
        # Compute residual at the current field value and form the new iterate g(x) = x - F(x)
        residual = objective_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
        x₁ = target_field[field_indices...] - residual
        # Apply (optionally relaxed) update
        target_field[field_indices(indices)...] = relaxed_update(solver.relax, x₁, x₀)
        iteration += 1
    end
    return x₁, iteration
end
