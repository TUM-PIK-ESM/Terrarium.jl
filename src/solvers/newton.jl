"""
    $TYPEDEF

Newton root-finder that always performs a *fixed* number of iterations, with the iteration count
carried in the type so that the iteration is a compile-time constant and the loop unrolls into
straight-line code.

Unlike [`RootSolver`](@ref) and [`FixedPointSolver`](@ref) there is no convergence test and hence no
data-dependent loop bound. That matters in two places:

- **Reactant**: a convergence-tested loop lowers to an `scf.while` with a dynamic trip count, which the
  StableHLO raise pass cannot lift. An unrolled loop raises cleanly.
- **Reverse-mode AD**: a fixed trip count avoids a dynamically-sized tape.

TODO: In the future, Reactant should just use the regular solvers. 
"""
struct NewtonSolver{NF, iterations} end

"""
    $TYPEDSIGNATURES

Construct a [`NewtonSolver`](@ref) performing exactly `iterations` Newton steps.
"""
function NewtonSolver(::Type{NF}; iterations::Int = 5) where {NF}
    iterations > 0 || throw(ArgumentError("iterations must be positive, got $iterations"))
    return NewtonSolver{NF, iterations}()
end

"""
    $TYPEDSIGNATURES

Number of Newton iterations performed by the given solver.
"""
@inline iterations(::NewtonSolver{NF, N}) where {NF, N} = N

@propagate_inbounds function solve!(
        out, indices::NTuple{M, Integer}, grid, fields,
        objective_func!::ObjectiveFunction{target},
        solver::NewtonSolver{NF, N},
        func_args...;
        func_kwargs...
    ) where {NF, M, N, target}
    target_field = getproperty(out, target)
    # `Field`s must always be indexed in 3D; see `field_indices`
    target_indices = field_indices(indices)
    # Forward-difference step size, scaled with the magnitude of the iterate
    h = cbrt(eps(NF))
    x = target_field[target_indices...]
    for _ in 1:N
        # The objective reads the target field, so it must be set before each evaluation
        target_field[target_indices...] = x
        residual = objective_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
        δ = h * (one(NF) + abs(x))
        target_field[target_indices...] = x + δ
        residual_perturbed = objective_func!(out, indices..., grid, fields, func_args...; func_kwargs...)
        ∂residual∂x = (residual_perturbed - residual) / δ
        # A vanishing derivative would send the step to ±Inf; hold the iterate instead. Note that
        # `ifelse` evaluates both branches, but the division itself is total (no throw path).
        step = ifelse(iszero(∂residual∂x), zero(NF), residual / ∂residual∂x)
        x = x - step
    end
    target_field[target_indices...] = x
    return x, N
end
