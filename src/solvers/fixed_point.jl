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
        step_func!::ObjectiveFunction{target},
        solver::FixedPointSolver{NF},
        func_args...;
        func_kwargs...
    ) where {NF, N, target}
    target_field = getproperty(out, target)
    x₀ = target_field[indices...]
    step_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
    x₁ = target_field[indices...]
    iteration = 1
    while abs(x₁ - x₀) > solver.tolerance && iteration <= solver.max_iterations
        x₀ = x₁
        # Compute objective function
        step_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
        # Read out result
        x₁ = target_field[indices...]
        # Apply update
        target_field[indices...] = relaxed_update(solver.relax, x₁, x₀)
        iteration += 1
    end
    return x₁, iteration
end
