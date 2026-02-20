using Terrarium
using Test

using Oceananigans.AbstractOperations: AbstractOperation

@testset "Root distribution" begin
    grid = ColumnGrid(CPU(), UniformSpacing(Δz = 0.1, N = 10))
    rootdist = StaticExponentialRootDistribution(Float64)
    vars = Terrarium.Variables(rootdist)
    state = initialize(vars, grid)
    # check that root fraction is defined as an operation
    @test isa(state.root_fraction, AbstractOperation)
    # check that the root fractions sum to one
    @test all(sum(state.root_fraction, dims = 3) .≈ one(eltype(grid)))
end
