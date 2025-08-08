using Terra
using Test

@testset "Autotrophic respiration" begin
    include("autotrophic_respiration_tests.jl")
end

@testset "Vegetation carbon dynamics" begin
    include("carbon_dynamics_tests.jl")
end

@testset "Vegetation dynamics" begin
    include("vegetation_dynamics_tests.jl")
end

@testset "Phenology" begin
    include("phenology_tests.jl")
end