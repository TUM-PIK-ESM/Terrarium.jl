using Terrarium
using Test

using Terrarium: tuplejoin, merge_duplicates, safediv, fastmap

@testset "Utilities" begin
    # tuplejoin
    @test tuplejoin() == ()
    @test tuplejoin((1,)) == (1,)
    @test tuplejoin((1,), (2,)) == (1,2)
    @test tuplejoin((1,),(1,)) == (1,1)
    @test tuplejoin((1,2),(3,4,5),(6,7,)) == (1,2,3,4,5,6,7)

    # merge_duplicates
    @test merge_duplicates(()) == ()
    @test merge_duplicates((1,)) == (1,)
    @test merge_duplicates((1,1,)) == (1,)
    @test merge_duplicates((1,2,2,3,4)) == (1,2,3,4)

    # safediv
    @test safediv(1, 2) ≈ 1/2
    @test safediv(1.0, 0.0) == 0
    @test safediv(0.0, 0.0) == 0
    @test safediv(0, 1) == 0

    # fastmap
    @test fastmap(x -> x + 1, (1,2,3)) == (2,3,4)
    @test fastmap(+, (1,2,3), (2,3,4)) == (3,5,7)
    @test fastmap(+, (a=1,b=2,c=3), (b=3,c=4,a=2)) == (a=3,b=5,c=7)
    @inferred fastmap(*, ("a", 3.0, 1), ("b", 2.0, missing))

    # Piecewise linear interpolation
    f = piecewise_linear(
        1.0u"m" => 1.0,
        0.0u"m" => -1.0,
        -2.0u"m" => -2.0
    )
    @test f(2.0) ≈ 1.0
    @test f(1.0) ≈ 1.0
    @test f(0.5) ≈ 0.0
    @test f(-1.0) ≈ -1.5
    @test f(-3.0) ≈ -2.0
end
