using DeltaLand
using DeltaLand: compute_γv, compute_ν_star
using Test

@testset "γv test" begin
    vegetation_dynamics = PALADYNVegetationDynamics()
    # For now, test that γv = γv_min 
    γv = compute_γv(vegetation_dynamics)
    @test γv == vegetation_dynamics.γv_min
end

@testset "ν_star" begin
    vegetation_dynamics = PALADYNVegetationDynamics()
    # Test ν < ν_seed (ν_star should be equal to ν_seed)
    ν = vegetation_dynamics.ν_seed / 2
    ν_star = compute_ν_star(vegetation_dynamics, ν)
    @test ν_star == vegetation_dynamics.ν_seed

    # Test ν > ν_seed (ν_star should be equal to ν)
    ν = vegetation_dynamics.ν_seed * 2
    ν_star = compute_ν_star(vegetation_dynamics, ν)
    @test ν_star == ν

    # Test ν == ν_seed (ν_star should be equal to ν_seed)
    ν = vegetation_dynamics.ν_seed
    ν_star = compute_ν_star(vegetation_dynamics, ν)
    @test ν_star == vegetation_dynamics.ν_seed
end



