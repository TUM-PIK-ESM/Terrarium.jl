using Terrarium
using Terrarium: compute_γv, compute_ν_star, compute_ν_tendency
using Test

@testset "γv test" begin
    veg_dynamics = PALADYNVegetationDynamics()
    # For now, test that γv = γv_min
    γv = compute_γv(veg_dynamics)
    @test γv == veg_dynamics.γv_min
end

@testset "ν_star test" begin
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

    # TODO ν should be between 0 and 1!

end
compute_ν_tendency
@testset "ν_tend test" begin
    veg_dynamics = PALADYNVegetationDynamics()
    vegcarbon_dynamics = PALADYNCarbonDynamics()
    # Test ν_tendency should be finite (could be positive or negative)
    LAI_b = (vegcarbon_dynamics.LAI_min + vegcarbon_dynamics.LAI_max) / 2.0
    C_veg = 0.5 # Mock value
    ν = 0.3 # Mock value
    ν_tendency = compute_ν_tendency(veg_dynamics, vegcarbon_dynamics, LAI_b, C_veg, ν)
    @test isfinite(ν_tendency)
end
