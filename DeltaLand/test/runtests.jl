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

@testset "Soil model and processes" begin
    include("soil/soil_model_tests.jl")
end