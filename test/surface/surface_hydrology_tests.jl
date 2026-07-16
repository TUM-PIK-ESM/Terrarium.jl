using Terrarium
using Test

@testset "Canopy hydrology" begin
    include("canopy_interception/canopy_interception_tests.jl")
end

@testset "Evapotranspiration" begin
    include("evapotranspiration/evapotranspiration_tests.jl")
end

@testset "Surface runoff" begin
    include("runoff/surface_runoff_tests.jl")
end
