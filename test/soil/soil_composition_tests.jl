using Terrarium
using Test

using Terrarium: porosity, saturation, liquid_fraction, organic_fraction, mineral_texture

@testset "SoilComposition" begin
    let por = 0.3,
            sat = 0.5,
            liq = 0.5,
            org = 0.5,
            texture = SoilTexture(sand = 0.5, clay = 0.5)
        soil = SoilComposition(; porosity = por, saturation = sat, liquid = liq, organic = org, texture)
        @test porosity(soil) == por
        @test saturation(soil) == sat
        @test liquid_fraction(soil) == liq
        @test organic_fraction(soil) == org
        @test mineral_texture(soil) == texture
    end

    # check that argument bounds are enforced
    @test_throws AssertionError SoilComposition(porosity = -1.0)
    @test_throws AssertionError SoilComposition(porosity = 2.0)
    @test_throws AssertionError SoilComposition(saturation = -1.0)
    @test_throws AssertionError SoilComposition(saturation = 2.0)
    @test_throws AssertionError SoilComposition(liquid = -1.0)
    @test_throws AssertionError SoilComposition(liquid = 2.0)
    @test_throws AssertionError SoilComposition(organic = -1.0)
    @test_throws AssertionError SoilComposition(organic = 2.0)
end

@testset "volumetric_fractions" begin
    let por = 0.3,
            sat = 0.5,
            liq = 0.5,
            org = 0.5,
            texture = SoilTexture(sand = 0.5, clay = 0.5)
        soil = SoilComposition(; porosity = por, saturation = sat, liquid = liq, organic = org, texture)
        fracs = volumetric_fractions(soil)
        @test fracs.water == por * sat * liq
        @test fracs.ice == por * sat * (1 - liq)
        @test fracs.air == por * (1 - sat)
        @test fracs.organic == (1 - por) * org
        @test fracs.mineral == (1 - por) * (1 - org)
    end
end
