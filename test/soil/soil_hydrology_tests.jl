using Terrarium
using Test

@testset "Hydraulic properties (prescribed)" begin
    # For prescribed hyraulic properties, just check that the returned values
    # match what was set.
    hydraulic_props = PrescribedHydraulics(
        cond_sat = 1e-6,
        porosity = 0.3,
        field_capacity = 0.1,
        wilting_point = 0.02,
    )
    @test saturated_hydraulic_conductivity(hydraulic_props) == hydraulic_props.cond_sat
    @test mineral_porosity(hydraulic_props) == hydraulic_props.porosity
    @test field_capacity(hydraulic_props) == hydraulic_props.field_capacity
    @test wilting_point(hydraulic_props) == hydraulic_props.wilting_point
end

@testset "Hydraulic properties (SURFEX)" begin
    # mineral porosity
    # TODO: should the hydraulic properties struct constructors also enforce parameter bounds?
    hydraulic_props = SURFEXHydraulics()
    por0 = mineral_porosity(hydraulic_props, SoilTexture(sand=0.0, silt=0.7, clay=0.3))
    @test por0 ≈ hydraulic_props.porosity
    for sand in 0.1:0.1:1.0
        silt = (1 - sand)*0.7
        clay = (1 - sand)*0.3
        por = mineral_porosity(hydraulic_props, SoilTexture(; sand, silt, clay))
        # test that increasing sand content decreases the mineral porosity;
        # we could also reproduce the calculation here, but that seems a bit redundant
        # and of course depends on the choice of parameters.
        @test 0 < por < por0
    end

    # check that wilting point is equal to zero when there is no clay
    wp0 = wilting_point(hydraulic_props, SoilTexture(sand=0.5, silt=0.5, clay=0.0))
    @test iszero(wp0)
    for clay in 0.1:0.1:1.0
        sand = (1-clay)*0.7
        silt = (1-clay)*0.3
        wp = wilting_point(hydraulic_props, SoilTexture(; sand, silt, clay))
        @test 0 < wp < 1
    end

    # check that field capacity is equal to zero when there is no clay
    fc0 = field_capacity(hydraulic_props, SoilTexture(sand=0.5, silt=0.5, clay=0.0))
    @test iszero(fc0)
    for clay in 0.1:0.1:1.0
        sand = (1-clay)*0.7
        silt = (1-clay)*0.3
        fc = field_capacity(hydraulic_props, SoilTexture(; sand, silt, clay))
        @test 0 < fc < 1
    end
end
