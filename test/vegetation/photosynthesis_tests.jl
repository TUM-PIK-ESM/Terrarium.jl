using Terra
using Terra: compute_kinetic_parameters, compute_Γ_star, compute_PAR, compute_APAR, compute_pres_i
using Terra: compute_t_stress, compute_c1_c2, compute_Vc_max, compute_JE_JC, compute_β
using Terra: compute_Rd, compute_Ag, compute_And, compute_photosynthesis
using Test


@testset "Kinetic parameters test" begin
    photo = LUEPhotosynthesis()
    # Test kinetic parameters should be positive
    T_air = 20.0 # °C
    τ, Kc, Ko, = compute_kinetic_parameters(photo, T_air)
    @test τ > 0
    @test Kc > 0
    @test Ko > 0

    # Test temperature dependence, kinetic parameters should increase/decrease with T_air
    T_air_warm = 30.0 # °C
    τ_warm, Kc_warm, Ko_warm = compute_kinetic_parameters(photo, T_air_warm)
    @test τ_warm < τ  
    @test Kc_warm > Kc  
    @test Ko_warm > Ko  
end

@testset "Γ_star test" begin
    photo = LUEPhotosynthesis()
    # Test Γ_star should be positive 
    τ = 3000.0 # Value from test values above
    pres_O2 = 20.9e3 # Pa
    Γ_star = compute_Γ_star(photo, τ, pres_O2)
    @test Γ_star > 0

    # Test τ dependence, Γ_star inv. proportional to τ
    τ_warm = 2000.0 # Value from test values above
    Γ_star_warm = compute_Γ_star(photo, τ_warm, pres_O2)
    @test Γ_star_warm > Γ_star
end

@testset "PAR test" begin
    photo = LUEPhotosynthesis()
    # Test PAR should be positive for a positive swdown
    swdown = 50.0 # W/m²
    PAR = compute_PAR(photo, swdown)
    @test PAR > 0

    # Test swdown = 0 (PAR should be 0)  
    swdown = 0.0 # W/m²
    PAR = compute_PAR(photo, swdown)
    @test PAR == 0

    # Test linearity, PAR should be linear in swdown
    swdown = 50.0 # W/m²
    PAR_50 = compute_PAR(photo, swdown)
    swdown = 100.0
    PAR_100 = compute_PAR(photo, swdown)
    @test PAR_100 ≈ 2.0 * PAR_50
end

@testset "APAR test" begin
    photo = LUEPhotosynthesis()
    swdown = 50.0 # W/m²
    # Test APAR should be positive for a positive LAI
    LAI = 5.0
    APAR = compute_APAR(photo, swdown, LAI)
    @test APAR > 0.0

    # Test LAI = 0 (APAR should be 0)
    LAI = 0.0
    APAR = compute_APAR(photo, swdown, LAI)
    @test APAR == 0

    # Test LAI = Inf (APAR should be = α_a*PAR)
    LAI = Inf
    APAR = compute_APAR(photo, swdown, LAI)
    @test APAR == photo.αa * compute_PAR(photo, swdown)
end

@testset "press_i test" begin
    photo = LUEPhotosynthesis()
    # Test λc = 0 (pi should be 0)
    λc = 0.0
    pres_a = 40.0 # Pa
    pres_i = compute_pres_i(photo, λc, pres_a)
    @test  pres_i == 0.0

    # Test λc = 1 (pi should be pa)
    λc = 1.0
    pres_i = compute_pres_i(photo, λc, pres_a)
    @test pres_i == pres_a

    # Test pi should be positive and less than pres_a for positive λc 
    λc = 0.5
    pres_a = 40.0 # Pa
    pres_i = compute_pres_i(photo, λc, pres_a)
    @test 0 < pres_i < pres_a
end

@testset "t_stress test" begin
    photo = LUEPhotosynthesis()
    # Test T_air < T_CO2_low (t_stress should be 0)
    T_air = photo.t_CO2_low * 2 # photo.t_CO2_low is negative
    t_stress = compute_t_stress(photo, T_air)
    @test t_stress == 0.0

    # Test T_air = T_CO2_low (t_stress should be 0)
    T_air = photo.t_CO2_low + eps()
    t_stress = compute_t_stress(photo, T_air)
    @test t_stress == 0.0

    # Test T_air > T_CO2_high (t_stress should be 0)
    T_air = photo.t_CO2_high * 2
    t_stress = compute_t_stress(photo, T_air)
    @test t_stress == 0.0

    # Test T_air = T_CO2_high (t_stress should be 0)
    T_air = photo.t_CO2_high
    t_stress = compute_t_stress(photo, T_air)
    @test t_stress == 0.0

    # Test T_CO2_low < T_air < T_CO2_high (t_stress should be between 0 and 1)
    T_air = (photo.t_CO2_low + photo.t_CO2_high) / 2
    t_stress = compute_t_stress(photo, T_air)
    @test 0.0 < t_stress < 1.0
end

@testset "c1 and c2 test" begin
    photo = LUEPhotosynthesis()
    T_air = 20.0 # °C
    Γ_star = 3.0 # Value from test values above
    Kc = 20.0 # Value from test values above
    Ko = 3.0e4 # Value from test values above
    pres_O2 = 20.9e3 # Pa

    # Test pi = Γ_star (c1 and c2 should be 0)
    pres_i = Γ_star
    c_1, c_2 = compute_c1_c2(photo, T_air, Γ_star, Kc, Ko, pres_i, pres_O2)
    @test c_1 == 0
    @test c_2 == 0

    # Test pi < Γ_star (c1 and c2 should be negative, c1 can be 0 if t_stress=0)
    pres_i = Γ_star / 2.0
    c_1, c_2 = compute_c1_c2(photo, T_air, Γ_star, Kc, Ko, pres_i, pres_O2)
    @test c_1 <= 0
    @test c_2 < 0

    # Test pi > Γ_star (c1 and c2 should be positive, c1 can be 0 if t_stress=0)
    pres_i = Γ_star * 2.0
    c_1, c_2 = compute_c1_c2(photo, T_air, Γ_star, Kc, Ko, pres_i, pres_O2)
    @test c_1 >= 0
    @test c_2 > 0
end

@testset "Vc_max test" begin
    photo = LUEPhotosynthesis()
    c_1 = 0.5 
    Kc = 20.0 # Value from test values above
    Ko = 3.0e4 # Value from test values above
    Γ_star = 3.0 # Value from test values above
    pres_i = 20.0 # Pa
    pres_O2 = 20.9e3 # Pa

    # Test APAR = 0 (Vc_max should be 0)
    APAR = 0.0
    Vc_max = compute_Vc_max(photo, c_1, APAR, Kc, Ko, Γ_star, pres_i, pres_O2)
    @test Vc_max == 0.0

    # Test Vc_max should be finite
    APAR = 5.0 # Value from test values above
    Vc_max = compute_Vc_max(photo, c_1, APAR, Kc, Ko, Γ_star, pres_i, pres_O2)
    @test isfinite(Vc_max)
end

@testset "JE and JC test" begin
    photo = LUEPhotosynthesis()
    # Test c1 = c2 = 0 (JE and JC should be 0) 
    c_1 = 0.0
    c_2 = 0.0
    APAR = 5.0 # Value from test values above
    Vc_max = 8.0 # Value from test values above
    JE, JC = compute_JE_JC(photo, c_1, c_2, APAR, Vc_max)
    @test JE == 0.0
    @test JC == 0.0

    # Test APAR = Vcmax = 0 (JE and JC should be 0)
    c_1 = 0.5
    c_2 = 0.5
    APAR = 0.0
    Vc_max = 0.0
    JE, JC = compute_JE_JC(photo, c_1, c_2, APAR, Vc_max)
    @test JE == 0.0
    @test JC == 0.0

    # Test JE, JC should be finite for the other cases
    c_1 = 0.5
    c_2 = 0.5
    APAR = 5.0 # Value from test values above
    Vc_max = 8.0 # Value from test values above
    JE, JC = compute_JE_JC(photo, c_1, c_2, APAR, Vc_max)
    @test isfinite(JE)
    @test isfinite(JC)
end

@testset "β test" begin
    photo = LUEPhotosynthesis()
    # For now, test β should be 1
    β = compute_β(photo)
    @test β == 1.0
end

@testset "Rd test" begin
    photo = LUEPhotosynthesis()
    # Test β=0 (Rd should be 0)
    β = 0.0
    Vc_max = 8.0 # Value from test values above
    Rd = compute_Rd(photo, Vc_max, β)
    @test Rd == 0.0

    # Test β=1 (Rd should be = α_C3*Vc_max)
    β = 1.0
    Vc_max = 8.0 # Value from test values above
    Rd = compute_Rd(photo, Vc_max, β)
    @test Rd == photo.α_C3 * Vc_max

    # Test Rd should be between 0 and α_C3*Vc_max for 0<β<1
    β = 0.5
    Vc_max = 8.0 # Value from test values above
    Rd = compute_Rd(photo, Vc_max, β)
    @test 0.0 < Rd < photo.α_C3 * Vc_max
end

@testset "Ag test" begin
    photo = LUEPhotosynthesis()
    c_1 = 0.5
    c_2 = 0.5
    APAR = 5.0 # Value from test values above
    Vc_max = 10.0 # Value from test values above

    # Test with β = 0 (Ag should be 0)
    β = 0.0
    Ag = compute_Ag(photo, c_1, c_2, APAR, Vc_max, β)
    @test Ag == 0

    # Test Ag should be finite for all other cases
    β = 0.5
    Ag = compute_Ag(photo, c_1, c_2, APAR, Vc_max, β)
    @test isfinite(Ag)   
end

@testset "And test" begin
    photo = LUEPhotosynthesis()
    c_1 = 0.5
    c_2 = 0.5
    β = 0.5
    APAR = 5.0 # Value from test values above
    Vc_max = 10.0 # Value from test values above
    Rd = 0.3 # Value from test values above
    And = compute_And(photo, c_1, c_2, APAR, Vc_max, β, Rd)
    @test isfinite(And)
end

@testset "Photosynthesis (GPP and Rd) test" begin
    photo = LUEPhotosynthesis()
    swdown = 200.0 # W/m²
    pres = 1.0e5 # Pa
    co2 = 400.0  # ppm
    λc = 0.5

    # Test T_air < -3 (GPP and Rd should be 0)
    T_air = -5.0 # °C
    LAI= 5.0
    GPP, Rd = compute_photosynthesis(photo, T_air, swdown, pres, co2, LAI, λc)
    @test GPP == 0.0
    @test Rd == 0.0

    # Test T_air > -3 and LAI=0 (GPP and Rd should be 0)
    T_air = 20.0 # °C
    LAI= 0.0
    GPP, Rd = compute_photosynthesis(photo, T_air, swdown, pres, co2, LAI, λc)
    @test GPP == 0.0
    @test Rd == 0.0

    # Test T_air > -3 and LAI > 0 (GPP and Rd should be positive)
    T_air = 20.0 # °C
    LAI = 5.0
    GPP, Rd = compute_photosynthesis(photo, T_air, swdown, pres, co2, LAI, λc)
    @test GPP > 0.0
    @test Rd > 0.0
end