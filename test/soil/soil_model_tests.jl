using Terrarium
using Test

@testset "Soil energy" begin
    include("soil_energy_tests.jl")
end

@testset "Soil hydrology" begin
    include("soil_hydrology_tests.jl")
end
