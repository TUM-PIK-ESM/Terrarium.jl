using Terrarium
using Test

@testset "SoilPorosity" begin
    # mineral porosity
    # TODO: should the hydraulic properties struct constructors also enforce parameter bounds?
    porosity = SoilPorositySURFEX(Float64)
    por0 = mineral_porosity(porosity, SoilTexture(sand = 0.0, silt = 0.7, clay = 0.3))
    @test por0 ≈ porosity.porosity_default
    for sand in 0.1:0.1:1.0
        silt = (1 - sand) * 0.7
        clay = (1 - sand) * 0.3
        por = mineral_porosity(porosity, SoilTexture(; sand, silt, clay))
        # test that increasing sand content decreases the mineral porosity;
        # we could also reproduce the calculation here, but that seems a bit redundant
        # and of course depends on the choice of parameters.
        @test 0 < por < por0
    end
end
