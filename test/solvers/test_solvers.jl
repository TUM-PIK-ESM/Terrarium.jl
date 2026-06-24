using Terrarium
using Test

using Terrarium: FixedPointSolver, RootSolver, RootSolvers, ObjectiveFunction, RelaxationFactor, solve!

# Build an objective function from a scalar map g(x). The solvers now expect the
# objective to *return* the fixed-point residual F(x) = x - g(x) evaluated at the
# current value of the target field, rather than writing g(x) into the field; the
# solver itself owns the field update. Both solver families solve x = g(x).
# The `grid` and `fields` arguments are unused for these scalar tests.
function scalar_objective(g, ::Val{target}) where {target}
    function residual(out, i, _grid, _fields)
        x = out.x[i]
        return x - g(x)
    end
    return ObjectiveFunction(residual, :x)
end

# Run a solver on g starting from x0 and return (root, iterations).
function run_solver(solver, g, x0)
    out = (; x = [x0])
    objective = scalar_objective(g, Val{:x}())
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
            x, iters = run_solver(solver, x -> x / 2 + 1, NF(0.0))
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
            x, iters = run_solver(solver, x -> (x + 2 / x) / 2, NF(1.0))
            @test isapprox(x, sqrt(NF(2.0)); atol = tol)
            @test iters <= solver.max_iterations
        end

        @testset "Iteration limit" begin
            # g(x) = 2x is non-contractive; the iteration cannot converge and must
            # exhaust the iteration budget.
            limited = FixedPointSolver(NF; tolerance = tol, relax = RelaxationFactor(factor = NF(1.0)), max_iterations = 5)
            _, iters = run_solver(limited, x -> 2 * x, NF(1.0))
            @test iters > limited.max_iterations
        end

        @testset "Type stability" begin
            # The full solve path must infer to a concrete (NF, Int) tuple with no
            # boxing (see the header note on avoiding `NF` capture in objectives).
            result = @inferred run_solver(solver, x -> x / 2 + 1, NF(0.0))
            @test result isa Tuple{NF, Int}
        end
    end

    @testset "RootSolver ($NF)" begin
        tolerance = RootSolvers.ResidualTolerance(tol)
        solver = RootSolver(NF; tolerance, max_iterations = 100)

        @testset "Linear map" begin
            # Same linear fixed point x = 2 as above
            x = run_solver(solver, x -> x / 2 + 1, NF(0.0))
            @test isapprox(x, NF(2.0); atol = tol)
        end

        @testset "Dottie number (cos)" begin
            x = run_solver(solver, cos, NF(1.0))
            @test isapprox(x, NF(0.7390851332151607); atol = tol)
        end

        @testset "Square root via x^2 - 2 = 0" begin
            # The fixed-point form g(x) = x^2 + x - 2 corresponds to the residual
            # F(x) = x - g(x) = 2 - x^2, whose positive root is √2.
            x = run_solver(solver, x -> x^2 + x - 2, NF(1.0))
            @test isapprox(x, sqrt(NF(2.0)); atol = tol)
        end

        @testset "Iteration limit" begin
            # A stiff, slowly-converging residual paired with a tight iteration budget:
            # the solver must stop after exactly `max_iterations` steps.
            limited_solver = RootSolver(NF; tolerance, max_iterations = 2, solution_type = RootSolvers.VerboseSolution())
            _, iters = run_solver(limited_solver, x -> (x - 10)^2 + log(cbrt(x + 2) + exp(x - 1)), NF(0.0))
            @test iters == limited_solver.max_iterations
        end

        @testset "Type stability" begin
            # With the default CompactSolution the solve path must infer to a
            # concrete scalar root of type NF with no boxing.
            x = @inferred run_solver(solver, x -> x / 2 + 1, NF(0.0))
            @test x isa NF
            @test isapprox(x, NF(2.0); atol = tol)
        end

        @testset "Verbose solution type" begin
            # VerboseSolution returns the iteration count alongside the root,
            # whereas the default CompactSolution returns only the scalar root.
            verbose_solver = RootSolver(
                NF;
                tolerance,
                solution_type = RootSolvers.VerboseSolution(),
                max_iterations = 100
            )
            x, iters = run_solver(verbose_solver, x -> x / 2 + 1, NF(0.0))
            @test isapprox(x, NF(2.0); atol = tol)
            @test iters isa Integer
            @test 0 < iters <= verbose_solver.max_iterations
        end
    end
end
