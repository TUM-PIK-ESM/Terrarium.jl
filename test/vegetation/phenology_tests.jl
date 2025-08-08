using Terra
using Terra: compute_f_deciduous, compute_phen
using Test

@testset "f_deciduous test" begin
    phenology = PALADYNPhenology{Float64}()
    # For now, test that f_deciduous = 0
    f_deciduous = compute_f_deciduous(phenology)
    @test f_deciduous == 0.0
end

@testset "phen test" begin
    phenology = PALADYNPhenology{Float64}()
    # For now, test that phen = 1
    phen = compute_phen(phenology)
    @test phen == 1.0
end
