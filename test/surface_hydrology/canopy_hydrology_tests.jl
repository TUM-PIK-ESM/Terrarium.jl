using Terrarium
using Terrarium:
    compute_canopy_interception,
    compute_canopy_saturation_fraction,
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

@testset "compute_w_can_tend" begin
    canopy_hydrology = PALADYNCanopyHydrology()
    constants = PhysicalConstants()
    # Test tendency is zero when there is no water present or incoming
    ∂w∂t, R_can = compute_w_can_tend(canopy_hydrology, constants, 0.0, 0.0, 0.0)
    @test iszero(∂w∂t)
    @test iszero(R_can)
    # Test tendency is negative when water is present but no evaporation or interception
    w_can = 0.1
    ∂w∂t, R_can = compute_w_can_tend(canopy_hydrology, constants, w_can, 0.0, 0.0)
    @test ∂w∂t == -R_can
    @test ∂w∂t < 0
    # Test interception and evaporation cancel 
    ∂w∂t, R_can = compute_w_can_tend(canopy_hydrology, constants, 0.0, 1e-6, 1e-6)
    @test iszero(∂w∂t)
    @test iszero(R_can)
    # Test positive with incoming interception
    ∂w∂t, R_can = compute_w_can_tend(canopy_hydrology, constants, 0.0, 1e-6, 1e-7)
    @test ∂w∂t > 0
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
