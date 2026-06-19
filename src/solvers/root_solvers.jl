import RootSolvers
import FiniteDifferences

"""
    $TYPEDEF

Wrapper for RootSolvers.jl root-finding methods.
"""
struct RootSolver{NF, M, S}
    "Numerical tolerance of the root finding iteration"
    tolerance::NF

    "Maximum number of iterations to run"
    max_iterations::Int

    "Solution type (e.g., VerboseSolution, CompactSolution)"
    solution_type::S

    function RootSolver(::Type{NF}, ::Type{M}, tolerance, max_iterations, solution_type) where {NF, M <: RootSolvers.RootSolvingMethod}
        return new{NF, M, typeof(solution_type)}(tolerance, max_iterations, solution_type)
    end
end

function RootSolver(
        ::Type{NF};
        method = RootSolvers.NewtonsMethod{NF},
        tolerance::NF = sqrt(eps(NF)),
        max_iterations::Int = 100,
        solution_type::S = RootSolvers.VerboseSolution()
    ) where {NF, S}
    return RootSolver(NF, method, tolerance, max_iterations, solution_type)
end

@propagate_inbounds function solve!(
        out, indices::NTuple{N, Integer}, grid, fields,
        step_func!::ObjectiveFunction{target},
        solver::RootSolver{NF, M, S},
        func_args...;
        func_kwargs...
    ) where {NF, N, M, S, target}
    target_field = getproperty(out, target)

    # Define the residual function F(x) = g(x) - x
    residual! = build_residual(out, indices, grid, fields, step_func!, target_field, solver, func_args...; func_kwargs...)

    # Solve for the root using RootSolvers.jl
    x₀ = target_field[indices...]
    method = M(x₀)
    tol = RootSolvers.ResidualTolerance(solver.tolerance)
    y₀ = residual!(x₀)
    result = RootSolvers.find_zero(residual!, method, solver.solution_type, tol, solver.max_iterations)

    # Extract the root and number of iterations from the result
    x_root = result.root
    iteration = S <: RootSolvers.CompactSolution ? nothing : result.iter_performed

    # Ensure the target field holds the final estimate
    target_field[indices...] = x_root

    return x_root, iteration
end


@propagate_inbounds function build_residual(
        out, indices, grid, fields,
        step_func!::ObjectiveFunction{<:Any, F, DF},
        target_field,
        solver::RootSolver{NF, RootSolvers.NewtonsMethod{NF}},
        func_args...;
        func_kwargs...
    ) where {NF, F, DF}

    dstep_func! = step_func!.dfunc

    function residual!(x)
        target_field[indices...] = x
        step_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
        gx = target_field[indices...]
        return gx - x
    end

    function dresidual!(x)
        target_field[indices...] = x
        dstep_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
        dgx = target_field[indices...]
        return dgx - 1
    end

    return x -> (residual!(x), dresidual!(x))
end

@propagate_inbounds function build_residual(
        out, indices, grid, fields,
        step_func!::ObjectiveFunction{<:Any, F, Nothing},
        target_field,
        solver::RootSolver{NF, RootSolvers.NewtonsMethod{NF}},
        func_args...;
        func_kwargs...
    ) where {NF, F}

    function residual!(x)
        target_field[indices...] = x
        step_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
        gx = target_field[indices...]
        return gx - x
    end

    fd = FiniteDifferences.central_fdm(2, 1)
    return x -> (residual!(x), fd(residual!, x)::NF)
end
