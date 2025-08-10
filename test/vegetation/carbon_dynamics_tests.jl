using Terra
using Terra: compute_λ_NPP, compute_LAI_b, compute_Λ_loc, compute_C_veg_tend
using Test

@testset "λ_NPP test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test LAI_b < LAI_min (λ_NPP should be 0)
    LAI_b = vegcarbon_dynamics.LAI_min / 2.0
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) 
    @test λ_NPP == 0.0

    # Test LAI_b == LAI_min (λ_NPP should be 0)
    LAI_b = vegcarbon_dynamics.LAI_min
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) 
    @test λ_NPP == 0.0

    # Test LAI_min < LAI_b < LAI_max (λ_NPP should be between 0 and 1)
    LAI_b = (vegcarbon_dynamics.LAI_min + vegcarbon_dynamics.LAI_max) / 2.0
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b)
    @test 0 < λ_NPP < 1

    # Test LAI_b == LAI_max (λ_NPP should be 1)
    LAI_b = vegcarbon_dynamics.LAI_max
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b)
    @test λ_NPP == 1.0

    # Test LAI_b > LAI_max (λ_NPP should be 1)
    LAI_b = vegcarbon_dynamics.LAI_max * 2.0
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b)
    @test λ_NPP == 1.0
end

@testset "LAI_b test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test LAI_b should be positive for a positive C_veg
    C_veg = 0.5 # Mock value
    LAI_b = compute_LAI_b(vegcarbon_dynamics, C_veg)
    @test isfinite(LAI_b) && LAI_b > 0.0

    # TODO C_veg, LAI_b should be always positive!
end

@testset "Λ_loc test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test Λ_loc should be positive for a positive LAI_b
    LAI_b = (vegcarbon_dynamics.LAI_min + vegcarbon_dynamics.LAI_max) / 2.0
    Λ_loc = compute_Λ_loc(vegcarbon_dynamics, LAI_b)
    @test isfinite(Λ_loc) && Λ_loc > 0.0

    # TODO C_veg, LAI_b should be always positive!
end

@testset "C_veg tendency test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test C_veg_tendency should be finite (could be positive or negative)
    LAI_b = (vegcarbon_dynamics.LAI_min + vegcarbon_dynamics.LAI_max) / 2.0
    NPP = 0.5 # Mock value
    C_veg_tendency = compute_C_veg_tend(vegcarbon_dynamics, LAI_b, NPP)
    @test isfinite(C_veg_tendency)
end
