using DeltaLand
using Test

@testset "Grids" begin
    include("grids.jl")
end

@testset "State variables" begin
    include("state_variables.jl")
end

@testset "Timestepping" begin
    include("timestepping/run_full_model.jl")
end