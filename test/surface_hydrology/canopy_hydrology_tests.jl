using Terrarium
using Terrarium:
    compute_canopy_interception,
    compute_canopy_saturation_fraction,
    compute_canopy_water_removal,
    compute_w_can_tend,
    compute_precip_ground
using Test

@testset "compute_canopy_interception" begin
    canopy_hydrology = PALADYNCanopyHydrology()
    # Test interception is zero when there is no rain
    I_can = compute_canopy_interception(canopy_hydrology, 0.0, 1.0, 0.5)
    @test iszero(I_can)
    # Also when LAI and SAI are zero
    I_can = compute_canopy_interception(canopy_hydrology, 1.0, 0.0, 0.0)
    @test iszero(I_can)
    # Test interception is between zero and the precipitation rate
    P = 1e-8
    I_can = compute_canopy_interception(canopy_hydrology, P, 1.0, 0.5)
    @test 0 < I_can < P
end

@testset "compute_canopy_saturation_fraction" begin
    canopy_hydrology = PALADYNCanopyHydrology()
    # Test saturated fraction is zero when there is no water
    f_can = compute_canopy_saturation_fraction(canopy_hydrology, 0.0, 1.0, 0.5)
    @test iszero(f_can)
    # Also when LAI and SAI are zero
    f_can = compute_canopy_saturation_fraction(canopy_hydrology, 1.0, 0.0, 0.0)
    @test iszero(f_can)
    # Test saturation is between 0 and 1
    w_can = 0.1
    f_can = compute_canopy_saturation_fraction(canopy_hydrology, w_can, 1.0, 0.5)
    @test 0 < f_can < 1
    # Check that f_can decreases when we increase the total LAI + SAI
    f_can2 = compute_canopy_saturation_fraction(canopy_hydrology, w_can, 2.0, 1.0)
    @test f_can2 < f_can
end

@testset "compute_canopy_water_removal" begin
    canopy_hydrology = PALADYNCanopyHydrology()
    constants = PhysicalConstants()
    # Test flux is zero when there is no stored water
    ∂w∂t = compute_canopy_water_removal(canopy_hydrology, constants, 0.0)
    @test iszero(∂w∂t)
    # Test that flux is still zero when there is negative water (mass balance violation)
    ∂w∂t = compute_canopy_water_removal(canopy_hydrology, constants, -1.0)
    @test iszero(∂w∂t)
    # Test flux is positive when there is water
    ∂w∂t = compute_canopy_water_removal(canopy_hydrology, constants, 1.0)
    @test ∂w∂t > 0
end

@testset "compute_w_can_tend" begin
    canopy_hydrology = PALADYNCanopyHydrology()
    constants = PhysicalConstants()
    # Test tendency is zero when all flux terms are zero
    ∂w∂t = compute_w_can_tend(canopy_hydrology, 0.0, 0.0, 0.0)
    @test iszero(∂w∂t)
    # Test tendency is negative when removal is positive
    ∂w∂t = compute_w_can_tend(canopy_hydrology, 0.0, 0.0, 1.0)
    @test ∂w∂t < 0
    # Test that interception and evaporation cancel 
    ∂w∂t = compute_w_can_tend(canopy_hydrology, 1e-6, 1e-6, 0.0)
    @test iszero(∂w∂t)
    # Test positive with incoming interception
    ∂w∂t = compute_w_can_tend(canopy_hydrology, 1e-6, 1e-7, 1e-7)
    @test ∂w∂t ≈ 1e-6 - 2e-7
end

@testset "compute_precip_ground" begin
    canopy_hydrology = PALADYNCanopyHydrology()
    precip_ground = compute_precip_ground(canopy_hydrology, 0, 0, 0)
    @test iszero(precip_ground)
    # Test calculation of precip_ground
    P = 1e-8
    I_can = P / 2
    R_can = 1e-6
    precip_ground = compute_precip_ground(canopy_hydrology, P, I_can, R_can)
    @test precip_ground == P - I_can + R_can
end
