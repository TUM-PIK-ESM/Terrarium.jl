using DeltaLand
using DeltaLand: compute_λ_NPP, compute_LAI_b, compute_Λ_loc
using Test

@testset "λ_NPP test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test LAI_b < LAI_min (should return 0)
    LAI_b = vegcarbon_dynamics.LAI_min/2.0
    @test λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) == 0.0

    # Test LAI_b == LAI_min (should return 0)
    LAI_b = vegcarbon_dynamics.LAI_min
    @test λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) == 0.0

    # Test LAI_min < LAI_b < LAI_max (should return a value between 0 and 1)
    LAI_b = (vegcarbon_dynamics.LAI_min + vegcarbon_dynamics.LAI_max) / 2
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b)
    @test 0 < λ_NPP < 1

    # Test LAI_b == LAI_max (should return 1)
    LAI_b = vegcarbon_dynamics.LAI_max
    @test λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) == 1.0

    # Test LAI_b > LAI_max (should return 1)
    LAI_b = vegcarbon_dynamics.LAI_max * 2.0
    @test λ_NPP = compute_λ_NPP(vegcarbon_dynamics, LAI_b) == 1.0
end

@testset "LAI_b test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Generate random positive value for C_veg
    C_veg = rand() + eps() # Ensure C_veg is positive
    # Test LAI_b should be positive
    LAI_b = compute_LAI_b(vegcarbon_dynamics, C_veg)
    @test LAI_b > 0.0 
end

@testset "Λ_loc test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Generate random positive value for LAI_b
    LAI_b = rand() + eps() # Ensure LAI_b is positive
    # Test Λ_loc should be positive
    Λ_loc = compute_Λ_loc(vegcarbon_dynamics, LAI_b)
    @test Λ_loc > 0.0
end

