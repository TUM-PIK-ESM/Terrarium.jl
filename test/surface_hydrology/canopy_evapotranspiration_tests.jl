using Terrarium
using Terrarium:
    compute_transpiration,
    compute_evaporation_ground,
    compute_evaporation_canopy
using Test

@testset "compute_transpiration" begin
    canopy_ET = PALADYNCanopyEvapotranspiration(Float64)
    # Test transpiration is zero when there is no VPD
    rₐ = 100
    gw_can = 0.1
    Δq = 0.0
    E_trp = compute_transpiration(canopy_ET, Δq, rₐ, gw_can)
    @test iszero(E_trp)
    # Test transpiration is positive and finite for nonzero VPD
    Δq = 0.01
    E_trp = compute_transpiration(canopy_ET, Δq, rₐ, gw_can)
    @test isfinite(E_trp)
    @test E_trp > 0
    # Test transpiration is still finite but small when resistances are large
    Δq = 0.01
    E_trp1 = compute_transpiration(canopy_ET, Δq, rₐ, 0.0)
    @test isfinite(E_trp1)
    @test 0 < E_trp1 < 1e-8
    @test E_trp1 < E_trp
    E_trp2 = compute_transpiration(canopy_ET, Δq, 1000.0, 1e-4)
    @test isfinite(E_trp2)
    @test 0 < E_trp2 < 1e-6
end

@testset "compute_evaporation_ground" begin
    canopy_ET = PALADYNCanopyEvapotranspiration(Float64)
    # Test evaporation is zero when there is no VPD
    rₐ = 50
    rₑ = 100
    Δq = 0.0
    β = 1.0
    E_gnd = compute_evaporation_ground(canopy_ET, Δq, β, rₐ, rₑ)
    @test iszero(E_gnd)
    # Test evaporation is positive when there is VPD
    rₐ = 50
    rₑ = 100
    Δq = 0.001
    β = 1.0
    E_gnd = compute_evaporation_ground(canopy_ET, Δq, β, rₐ, rₑ)
    @test E_gnd > 0
    # Test evaporation is smaller when resistance is higher
    rₐ = 100
    rₑ = 100
    Δq = 0.001
    β = 1.0
    E_gnd2 = compute_evaporation_ground(canopy_ET, Δq, β, rₐ, rₑ)
    @test 0 < E_gnd2 < E_gnd
    # Test evaporation is smaller when limiting factor is < 1
    rₐ = 50
    rₑ = 100
    Δq = 0.001
    β = 0.5
    E_gnd2 = compute_evaporation_ground(canopy_ET, Δq, β, rₐ, rₑ)
    @test 0 < E_gnd2 < E_gnd
    @test E_gnd2 ≈ E_gnd / 2
    # Test evaporation is larger with larger VPD
    rₐ = 50
    rₑ = 100
    Δq = 0.01
    β = 1.0
    E_gnd3 = compute_evaporation_ground(canopy_ET, Δq, β, rₐ, rₑ)
    @test 0 < E_gnd2 < E_gnd < E_gnd3
end

@testset "compute_evaporation_canopy" begin
    canopy_ET = PALADYNCanopyEvapotranspiration(Float64)
    # Test evaporation is zero when there is no VPD
    rₐ = 50
    Δq = 0.0
    f_can = 1.0
    E_can = compute_evaporation_canopy(canopy_ET, Δq, f_can, rₐ)
    @test iszero(E_can)
    # Test evaporation is positive when there is VPD
    rₐ = 50
    Δq = 0.001
    f_can = 1.0
    E_can = compute_evaporation_canopy(canopy_ET, Δq, f_can, rₐ)
    @test E_can > 0
    # Test evaporation is smaller when resistance is higher
    rₐ = 100
    rₑ = 100
    Δq = 0.001
    f_can = 1.0
    E_can2 = compute_evaporation_canopy(canopy_ET, Δq, f_can, rₐ)
    @test 0 < E_can2 < E_can
    # Test evaporation is smaller when limiting factor is < 1
    rₐ = 50
    rₑ = 100
    Δq = 0.001
    f_can = 0.5
    E_can2 = compute_evaporation_canopy(canopy_ET, Δq, f_can, rₐ)
    @test 0 < E_can2 < E_can
    @test E_can2 ≈ E_can / 2
    # Test evaporation is larger with larger VPD
    rₐ = 50
    rₑ = 100
    Δq = 0.01
    f_can = 1.0
    E_can3 = compute_evaporation_canopy(canopy_ET, Δq, f_can, rₐ)
    @test 0 < E_can2 < E_can < E_can3
end
