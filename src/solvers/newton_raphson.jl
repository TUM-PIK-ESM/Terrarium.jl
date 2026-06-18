struct NewtonRaphsonSolver{NF, R}
    "Numerical tolerance of the Newton-Raphson iteration"
    tolerance::NF

    "Finite difference step size used to approximate the derivative"
    fd_step::NF

    "Relaxation scheme"
    relax::R

    "Maximum number of iterations to run"
    max_iterations::Int
end

function NewtonRaphsonSolver(
        ::Type{NF};
        tolerance::NF = sqrt(eps(NF)),
        fd_step::NF = cbrt(eps(NF)),
        relax::R = nothing, # default to no relaxation
        max_iterations::Int = 100
    ) where {NF, R}
    return NewtonRaphsonSolver{NF, R}(tolerance, fd_step, relax, max_iterations)
end

# Evaluate the fixed-point residual F(x) = g(x) - x, where g(x) is the value of
# the target field produced by `step_func!` when the field is initialized to `x`.
@propagate_inbounds function _residual!(
        out, indices, grid, fields, step_func!, target_field, x,
        func_args...; func_kwargs...
    )
    target_field[indices...] = x
    step_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
    gx = target_field[indices...]
    return gx - x
end

@propagate_inbounds function solve!(
        out, indices::NTuple{N, Integer}, grid, fields,
        step_func!::ObjectiveFunction{target},
        solver::NewtonRaphsonSolver{NF},
        func_args...;
        func_kwargs...
    ) where {NF, N, target}
    target_field = getproperty(out, target)
    x = target_field[indices...]
    # Evaluate residual at the initial guess
    F = _residual!(out, indices, grid, fields, step_func!, target_field, x, func_args...; func_kwargs...)
    iteration = 1
    while abs(F) > solver.tolerance && iteration <= solver.max_iterations
        # Approximate the derivative F'(x) with a forward finite difference
        h = solver.fd_step * (1 + abs(x))
        F_h = _residual!(out, indices, grid, fields, step_func!, target_field, x + h, func_args...; func_kwargs...)
        dF = (F_h - F) / h
        # Newton update; fall back to the perturbed point if the derivative vanishes
        x_new = iszero(dF) ? x + h : x - F / dF
        # Apply (optionally relaxed) update
        x = relaxed_update(solver.relax, x_new, x)
        # Re-evaluate residual at the updated estimate
        F = _residual!(out, indices, grid, fields, step_func!, target_field, x, func_args...; func_kwargs...)
        iteration += 1
    end
    # Ensure the target field holds the final estimate
    target_field[indices...] = x
    return x, iteration
end
