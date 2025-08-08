using Terra
using Terra: compute_f_temp, compute_resp10, compute_Rm, compute_Rg, compute_Ra, compute_NPP
using Test

@testset "f_temp test" begin
    auto_respiration = PALADYNAutotrophicRespiration()
    # For now, test that f_temp_soil= 0 and f_temp_air > 0
    T_air = 20.0
    f_temp_air, f_temp_soil = compute_f_temp(auto_respiration, T_air)
    @test f_temp_air > 0.0
    @test f_temp_soil == 0.0
end

@testset "resp10 test" begin
    auto_respiration = PALADYNAutotrophicRespiration()
    # For now, test that resp10 = 0.066
    resp10 = compute_resp10(auto_respiration)
    @test resp10 == 0.066
end

@testset "Rm test" begin
    auto_respiration = PALADYNAutotrophicRespiration()
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test Rm should be positive for a positive C_veg and Rd
    T_air = 20.0
    Rd = 500.0
    phen = 1.0
    C_veg = 0.5
    Rm = compute_Rm(auto_respiration, vegcarbon_dynamics, T_air, Rd, phen, C_veg)
    @test Rm > 0.0
end

@testset "Rg test" begin
    auto_respiration = PALADYNAutotrophicRespiration()
    # Test Rg should be finite
    GPP = 0.5
    Rm = 0.2
    Rg = compute_Rg(auto_respiration, GPP, Rm)
    @test isfinite(Rg)
end

@testset "Ra test" begin
    auto_respiration = PALADYNAutotrophicRespiration()
    # Test Ra should be the sum of Rm and Rg
    Rm = 0.2
    Rg = 0.1
    Ra = compute_Ra(auto_respiration, Rm, Rg)
    @test Ra == Rm + Rg
end

@testset "NPP test" begin
    auto_respiration = PALADYNAutotrophicRespiration()
    # Test NPP should be the difference between GPP and Ra
    GPP = 0.5
    Ra = 0.3
    NPP = compute_NPP(auto_respiration, GPP, Ra)
    @test NPP == GPP - Ra
end