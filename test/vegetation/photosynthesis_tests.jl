using Terra
using Terra: compute_t_stress, compute_kinetic_parameters, compute_PAR, compute_APAR, compute_pi
using Test

@testset "t_stress test" begin
    photosynthesis = LUEPhotosynthesis()
    # Test T_air < T_CO2_low (t_stress should be 0)
    T_air = photosynthesis.t_CO2_low * 2 # photosynthesis.t_CO2_low is negative
    t_stress = compute_t_stress(photosynthesis, T_air)
    @test t_stress == 0.0

    # Test T_air = T_CO2_low (t_stress should be 0)
    T_air = photosynthesis.t_CO2_low + eps()
    t_stress = compute_t_stress(photosynthesis, T_air)
    @test t_stress == 0.0

    # Test T_air > T_CO2_high (t_stress should be 0)
    T_air = photosynthesis.t_CO2_high * 2
    t_stress = compute_t_stress(photosynthesis, T_air)
    @test t_stress == 0.0

    # Test T_air = T_CO2_high (t_stress should be 0)
    T_air = photosynthesis.t_CO2_high
    t_stress = compute_t_stress(photosynthesis, T_air)
    @test t_stress == 0.0

    # Test T_CO2_low < T_air < T_CO2_high (t_stress should be between 0 and 1)
    T_air = (photosynthesis.t_CO2_low + photosynthesis.t_CO2_high) / 2
    t_stress = compute_t_stress(photosynthesis, T_air)
    @test 0.0 < t_stress < 1.0
end

@testset "Kinetic parameters test" begin
    photosynthesis = LUEPhotosynthesis()
    # Test kinetic parameters should be positive
    T_air = 20.0
    pres = 1.0e5
    τ, Kc, Ko, Γ_star = compute_kinetic_parameters(photosynthesis, T_air, pres)
    @test τ > 0
    @test Kc > 0
    @test Ko > 0
    @test Γ_star > 0
end

@testset "PAR test" begin
    photosynthesis = LUEPhotosynthesis()
    # Test PAR should be positive for a positive swdown
    swdown = 50.0
    PAR = compute_PAR(photosynthesis, swdown)
    @test PAR > 0

    # Test PAR should be 0 for swdown = 0
    swdown = 0.0
    PAR = compute_PAR(photosynthesis, swdown)
    @test PAR == 0
end

@testset "APAR test" begin
    photosynthesis = LUEPhotosynthesis()
    # Test APAR should be positive and less than PAR for a fixed LAI
    LAI = 5.0
    PAR = 30.0
    APAR = compute_APAR(photosynthesis, PAR, LAI)
    @test 0 < APAR < PAR

    # Test LAI = 0 (APAR should be 0)
    LAI = 0.0
    PAR = 30.0
    APAR = compute_APAR(photosynthesis, PAR, LAI)
    @test APAR == 0
end

@testset "pi test" begin
    photosynthesis = LUEPhotosynthesis()
    # Test for positive λc pi should be positive and less than pa
    λc = 0.5
    pa = 1.0e5
    pi = compute_pi(photosynthesis, λc, pa)
    @test 0 < pi < pa

    #TODO test for λc < 0, does this makes sense?
end