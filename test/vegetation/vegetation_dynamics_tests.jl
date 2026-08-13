using Terrarium
using Terrarium: compute_γv, compute_ν_star, compute_ν_tendency
using Test

@testset "compute_γv test" begin
    veg_dynamics = PALADYNVegetationDynamics()
    # For now, test that γv = γv_min
    γv = compute_γv(veg_dynamics)
    @test γv == veg_dynamics.γv_min / (365.25 * 24 * 3600) # in seconds r^-1
end

@testset "compute_ν_star test" begin
    veg_dynamics = PALADYNVegetationDynamics()
    # Test ν < ν_seed (ν_star should be equal to ν_seed)
    ν = veg_dynamics.ν_seed / 2
    ν_star = compute_ν_star(veg_dynamics, ν)
    @test ν_star == veg_dynamics.ν_seed

    # Test ν > ν_seed (ν_star should be equal to ν)
    ν = veg_dynamics.ν_seed * 2
    ν_star = compute_ν_star(veg_dynamics, ν)
    @test ν_star == ν

    # Test ν == ν_seed (ν_star should be equal to ν_seed)
    ν = veg_dynamics.ν_seed
    ν_star = compute_ν_star(veg_dynamics, ν)
    @test ν_star == veg_dynamics.ν_seed
end

@testset "compute_ν_tendency test" begin
    veg_dynamics = PALADYNVegetationDynamics()
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    traits = PlantTraits()
    # Test ν_tendency should be finite (could be positive or negative)
    LAI_b = (traits.minimum_leaf_area_index + traits.maximum_leaf_area_index) / 2.0
    # Mock values for state variables
    C_veg = 0.5
    ν = 0.3
    NPP = 1.0e-3
    ν_tendency = compute_ν_tendency(veg_dynamics, vegcarbon_dynamics, traits, LAI_b, C_veg, NPP, ν)
    @test isfinite(ν_tendency)
end
