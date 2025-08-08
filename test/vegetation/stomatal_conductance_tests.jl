using Terra
using Terra: compute_λc
using Test

@testset "λc test" begin
    stom_conductance = MedlynStomatalConductance()
    # Test vpd near zero (λc should be close to 1)
    vpd = eps() 
    λc = compute_λc(stom_conductance, vpd)
    @test λc ≈ 1.0

    # Test vpd near inf (λc should be close to 0)
    # TODO here is it ok to use Inf to test the limit case? or a very large number?
    vpd = Inf
    λc = compute_λc(stom_conductance, vpd)
    @test λc ≈ -0.6

    # Test λc should be between -0.6 and 1 for a positive vpd
    vpd = 500
    λc = compute_λc(stom_conductance, vpd)
    @test -0.6 < λc < 1.0
end
   
