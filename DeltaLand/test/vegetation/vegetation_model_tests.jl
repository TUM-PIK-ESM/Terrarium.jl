using DeltaLand
using Test

@testset "Vegetation carbon dynamics" begin
    include("carbon_dynamics_tests.jl")
end

@testset "Vegetation dynamics" begin
    include("vegetation_dynamics_tests.jl")
end

@testset "Phenology" begin
    include("phenology_tests.jl")
end