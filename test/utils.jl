using Terrarium
using Test
using Unitful

using Terrarium: tuplejoin, merge_duplicates, merge_recursive, safediv, fastmap, piecewise_linear

@testset "Utilities" begin
    # tuplejoin
    @test tuplejoin() == ()
    @test tuplejoin((1,)) == (1,)
    @test tuplejoin((1,), (2,)) == (1, 2)
    @test tuplejoin((1,), (1,)) == (1, 1)
    @test tuplejoin((1, 2), (3, 4, 5), (6, 7)) == (1, 2, 3, 4, 5, 6, 7)

    # merge_duplicates
    @test merge_duplicates(()) == ()
    @test merge_duplicates((1,)) == (1,)
    @test merge_duplicates((1, 1)) == (1,)
    @test merge_duplicates((1, 2, 2, 3, 4)) == (1, 2, 3, 4)

    # merge recursive
    @test merge_recursive((;), (;)) == merge((;), (;)) == (;)
    @test merge_recursive((;), (a = 1,)) == merge((;), (a = 1,)) == (a = 1,)
    @test merge_recursive((a = 1,), (;)) == merge((a = 1,), (;)) == (a = 1,)
    @test merge_recursive((a = 1,), (b = 2,)) == merge((a = 1,), (b = 2,)) == (a = 1, b = 2)
    @test merge_recursive((a = 1,), (a = 2,)) == merge((a = 1,), (a = 2,)) == (a = 2,)
    @test merge_recursive((a = 1, b = 2), (a = 2, c = 3)) == merge((a = 1, b = 2), (a = 2, c = 3)) == (a = 2, b = 2, c = 3)
    @test merge_recursive((a = 1, b = (x = 2, y = 3)), (a = 2, b = (x = 3, z = 4))) == (a = 2, b = (x = 3, y = 3, z = 4))
    @test merge_recursive(
        (a = 1, b = (x = (u = -1, v = -2), y = 3)),
        (a = 2, b = (x = (v = -3, w = 0), z = 4))
    ) == (a = 2, b = (x = (u = -1, v = -3, w = 0), y = 3, z = 4))
    @test merge_recursive((a = 1,), nothing) == (a = 1,)
    @test merge_recursive(nothing, (a = 1,)) == (a = 1,)
    @test_throws MethodError merge_recursive(nothing, nothing)

    # safediv
    @test safediv(1, 2) ≈ 1 / 2
    @test safediv(1.0, 0.0) == 0
    @test safediv(0.0, 0.0) == 0
    @test safediv(0, 1) == 0

    # fastmap
    @test fastmap(x -> x + 1, (1, 2, 3)) == (2, 3, 4)
    @test fastmap(+, (1, 2, 3), (2, 3, 4)) == (3, 5, 7)
    @test fastmap(+, (a = 1, b = 2, c = 3), (b = 3, c = 4, a = 2)) == (a = 3, b = 5, c = 7)
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
