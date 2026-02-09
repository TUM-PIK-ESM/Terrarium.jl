using Terrarium
using Test

@testset "Canopy hydrology" begin
    include("canopy_interception_tests.jl")
end

@testset "Canopy ET" begin
    include("canopy_evapotranspiration_tests.jl")
end

@testset "Surface runoff" begin
    include("surface_runoff_tests.jl")
end
