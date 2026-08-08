using Terrarium
using Terrarium: compute_λ_NPP, compute_balanced_leaf_area_index, compute_Λ_loc, compute_C_veg_tend
using Test

@testset "λ_NPP test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    traits = PlantTraits()
    # Test LAI_b < minimum_leaf_area_index (λ_NPP should be 0)
    LAI_b = traits.minimum_leaf_area_index / 2.0
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, traits, LAI_b)
    @test λ_NPP == 0.0

    # Test LAI_b == minimum_leaf_area_index (λ_NPP should be 0)
    LAI_b = traits.minimum_leaf_area_index
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, traits, LAI_b)
    @test λ_NPP == 0.0

    # Test minimum_leaf_area_index < LAI_b < maximum_leaf_area_index (λ_NPP should be between 0 and 1)
    LAI_b = (traits.minimum_leaf_area_index + traits.maximum_leaf_area_index) / 2.0
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, traits, LAI_b)
    @test 0 < λ_NPP < 1

    # Test LAI_b == maximum_leaf_area_index (λ_NPP should be 1)
    LAI_b = traits.maximum_leaf_area_index
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, traits, LAI_b)
    @test λ_NPP == 1.0

    # Test LAI_b > maximum_leaf_area_index (λ_NPP should be 1)
    LAI_b = traits.maximum_leaf_area_index * 2.0
    λ_NPP = compute_λ_NPP(vegcarbon_dynamics, traits, LAI_b)
    @test λ_NPP == 1.0
end

@testset "LAI_b test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    traits = PlantTraits()
    # Test LAI_b should be positive for a positive C_veg
    C_veg = 0.5 # Mock value
    LAI_b = compute_balanced_leaf_area_index(vegcarbon_dynamics, traits, C_veg)
    @test isfinite(LAI_b) && LAI_b > 0.0

    # TODO C_veg, LAI_b should be always positive!
end

@testset "Λ_loc test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    traits = PlantTraits()
    # Test Λ_loc should be positive for a positive LAI_b
    LAI_b = (traits.minimum_leaf_area_index + traits.maximum_leaf_area_index) / 2.0
    Λ_loc = compute_Λ_loc(vegcarbon_dynamics, traits, LAI_b)
    @test isfinite(Λ_loc) && Λ_loc > 0.0

    # TODO C_veg, LAI_b should be always positive!
end

@testset "C_veg tendency test" begin
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    traits = PlantTraits()
    # Test C_veg_tendency should be finite (could be positive or negative)
    LAI_b = (traits.minimum_leaf_area_index + traits.maximum_leaf_area_index) / 2.0
    NPP = 0.5 # Mock value
    C_veg_tendency = compute_C_veg_tend(vegcarbon_dynamics, traits, LAI_b, NPP)
    @test isfinite(C_veg_tendency)
end
