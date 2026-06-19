struct NewtonRaphsonSolver{NF, R}
    "Numerical tolerance of the Newton-Raphson iteration"
    tolerance::NF

    "Finite difference step size used to approximate the derivative"
    epsilon::NF

    "Relaxation scheme"
    relax::R

    "Maximum number of iterations to run"
    max_iterations::Int
end

function NewtonRaphsonSolver(
        ::Type{NF};
        tolerance::NF = sqrt(eps(NF)),
        epsilon::NF = cbrt(eps(NF)),
        relax::R = nothing, # default to no relaxation
        max_iterations::Int = 100
    ) where {NF, R}
    return NewtonRaphsonSolver{NF, R}(tolerance, epsilon, relax, max_iterations)
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
    # Build residual function
    residual! = build_residual(out, indices, grid, fields, step_func!, target_field, solver, func_args...; func_kwargs...)
    # Evaluate residual at the initial guess
    F = residual!(x)
    iteration = 1
    while abs(F) > solver.tolerance && iteration <= solver.max_iterations
        # Finite-difference step size scaled to the current iterate
        h = solver.epsilon * (1 + abs(x))
        # Approximate the derivative F'(x) with a forward finite difference
        F_h = residual!(x + h)
        dF = (F_h - F) / h
        # Newton update; fall back to the perturbed point if the derivative vanishes
        x_new = iszero(dF) ? x + h : x - F / dF
        # Apply (optionally relaxed) update
        x = relaxed_update(solver.relax, x_new, x)
        # Re-evaluate the residual at the new iterate for the convergence test and next step
        F = residual!(x)
        iteration += 1
    end
    # Ensure the target field holds the final estimate
    target_field[indices...] = x
    return x, iteration
end
