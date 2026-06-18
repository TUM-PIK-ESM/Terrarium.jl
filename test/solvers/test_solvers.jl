using Terrarium
using Test

using Terrarium: FixedPointSolver, NewtonRaphsonSolver, ObjectiveFunction, RelaxationFactor, solve!

# Build an objective function from a scalar map g(x). Both solvers treat `step_func!`
# as writing g(x) into the target field, i.e. they solve the fixed-point equation
# x = g(x). The `grid` and `fields` arguments are unused for these scalar tests.
function scalar_objective(g)
    func = function (out, i, _grid, _fields)
        out.x[i] = g(out.x[i])
        return nothing
    end
    return ObjectiveFunction(func, :x)
end

# Run a solver on g starting from x0 and return (root, iterations).
function run_solver(solver, g, x0)
    out = (; x = [x0])
    objective = scalar_objective(g)
    return solve!(out, (1,), nothing, nothing, objective, solver)
end

for NF in (Float32, Float64)
    tol = sqrt(eps(NF))

    @testset "FixedPointSolver ($NF)" begin
        # Use an undamped (factor = 1, or relax = nothing) update so the iteration matches the plain
        # fixed-point map for these contractive examples.
        solver = FixedPointSolver(NF; tolerance = tol, relax = RelaxationFactor(factor = NF(1.0)), max_iterations = 1000)

        @testset "Linear map" begin
            # g(x) = x/2 + 1 has the unique fixed point x = 2
            x, iters = run_solver(solver, x -> x / NF(2) + NF(1), NF(0.0))
            @test isapprox(x, NF(2.0); atol = tol)
            @test iters <= solver.max_iterations
        end

        @testset "Dottie number (cos)" begin
            # x = cos(x) ⇒ Dottie number ≈ 0.7390851332151607
            x, iters = run_solver(solver, cos, NF(1.0))
            @test isapprox(x, NF(0.7390851332151607); atol = tol)
            @test iters <= solver.max_iterations
        end

        @testset "Square root (Heron)" begin
            # g(x) = (x + 2/x)/2 has fixed point √2
            x, iters = run_solver(solver, x -> (x + NF(2) / x) / NF(2), NF(1.0))
            @test isapprox(x, sqrt(NF(2.0)); atol = tol)
            @test iters <= solver.max_iterations
        end

        @testset "Iteration limit" begin
            # g(x) = 2x is non-contractive; the iteration cannot converge and must
            # exhaust the iteration budget.
            limited = FixedPointSolver(NF; tolerance = tol, relax = RelaxationFactor(factor = NF(1.0)), max_iterations = 5)
            _, iters = run_solver(limited, x -> NF(2) * x, NF(1.0))
            @test iters > limited.max_iterations
        end
    end

    @testset "NewtonRaphsonSolver ($NF)" begin
        solver = NewtonRaphsonSolver(NF; tolerance = tol, max_iterations = 100)

        @testset "Linear map" begin
            # Same linear fixed point x = 2 as above
            x, iters = run_solver(solver, x -> x / NF(2) + NF(1), NF(0.0))
            @test isapprox(x, NF(2.0); atol = tol)
            @test iters <= solver.max_iterations
        end

        @testset "Dottie number (cos)" begin
            x, iters = run_solver(solver, cos, NF(1.0))
            @test isapprox(x, NF(0.7390851332151607); atol = tol)
            @test iters <= solver.max_iterations
        end

        @testset "Square root via x^2 - 2 = 0" begin
            # The fixed-point form g(x) = x^2 + x - 2 corresponds to the residual
            # F(x) = g(x) - x = x^2 - 2, whose positive root is √2. The plain
            # fixed-point iteration on this g diverges, but Newton converges.
            x, iters = run_solver(solver, x -> x^2 + x - NF(2), NF(1.0))
            @test isapprox(x, sqrt(NF(2.0)); atol = tol)
            @test iters <= solver.max_iterations
        end

        @testset "Fast convergence" begin
            # Newton's quadratic convergence should reach the root in only a few
            # iterations for a smooth, well-conditioned problem.
            _, iters = run_solver(solver, x -> x^2 + x - NF(2), NF(1.0))
            @test iters <= 10
        end

        @testset "Iteration limit" begin
            # F(x) = g(x) - x = 1 has no root, so the derivative vanishes and the
            # solver can never converge; it must exhaust the iteration budget.
            limited = NewtonRaphsonSolver(NF; tolerance = tol, max_iterations = 5)
            _, iters = run_solver(limited, x -> x + NF(1), NF(0.0))
            @test iters > limited.max_iterations
        end
    end
end
