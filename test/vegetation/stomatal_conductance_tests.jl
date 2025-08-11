using Terra
using Terra: compute_λc
using Test

@testset "λc test" begin
    stomcond = MedlynStomatalConductance()
    # Test vpd near zero (λc should be close to 1)
    vpd = eps() 
    λc = compute_λc(stomcond, vpd)
    @test λc ≈ 1.0

    # Test λc should be between 0 and 1 for a realistic positive vpd
    # TODO for very high VPD λc can become negative
    vpd = 1000 # Pa
    λc = compute_λc(stomcond, vpd)
    @test 0.0 < λc < 1.0
end
   
