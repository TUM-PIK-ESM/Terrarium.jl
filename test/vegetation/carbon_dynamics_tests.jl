using Terra
using Terra: compute_λ_NPP, compute_LAI_b, compute_Λ_loc
using Test

@testset "λ_NPP test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test LAI_b < LAI_min (λ_NPP should be 0)
    LAI_b = vegcarbon_dynamics.LAI_min/2.0
    @test λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) == 0.0

    # Test LAI_b == LAI_min (λ_NPP should be 0)
    LAI_b = vegcarbon_dynamics.LAI_min
    @test λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) == 0.0

    # Test LAI_min < LAI_b < LAI_max (λ_NPP should be between 0 and 1)
    LAI_b = (vegcarbon_dynamics.LAI_min + vegcarbon_dynamics.LAI_max) / 2
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b)
    @test 0 < λ_NPP < 1

    # Test LAI_b == LAI_max (λ_NPP should be 1)
    LAI_b = vegcarbon_dynamics.LAI_max
    @test λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) == 1.0

    # Test LAI_b > LAI_max (λ_NPP should be 1)
    LAI_b = vegcarbon_dynamics.LAI_max * 2.0
    @test λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) == 1.0
end

@testset "LAI_b test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test LAI_b should be positive for a positive C_veg
    C_veg = 0.5
    LAI_b = compute_LAI_b(vegcarbon_dynamics, C_veg)
    @test LAI_b > 0.0 
end

@testset "Λ_loc test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test Λ_loc should be positive for a positive LAI_b
    LAI_b = 0.5
    Λ_loc = compute_Λ_loc(vegcarbon_dynamics, LAI_b)
    @test Λ_loc > 0.0
end

