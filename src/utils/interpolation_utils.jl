"""
    $SIGNATURES

Return a function `f(z)` that linearly interpolates between the given `knots`.
"""
function piecewise_linear(knots::Pair{<:LengthQuantity}...; extrapolation = Interpolations.Flat())
    # extract coordinates and strip units
    zs = collect(map(ustrip ∘ first, knots))
    ys = collect(map(last, knots))
    @assert issorted(zs, rev = true) "depths must be sorted in descending order"
    interp = Interpolations.interpolate((reverse(zs),), reverse(ys), Interpolations.Gridded(Interpolations.Linear()))
    return Interpolations.extrapolate(interp, extrapolation)
end
