"""
    $SIGNATURES

Convert `Δt`s of type `Period` to a numeric value in seconds. Return `Δt` if already a number.
"""
convert_dt(Δt::Number) = Δt
convert_dt(Δt::Period) = Second(Δt).value

"""
    $SIGNATURES

Evaluates `x / (y + eps(NF))` if and only if `y != zero(y)`; returns `Inf` otherwise.
"""
safediv(x::NF, y::NF) where {NF} = ifelse(iszero(y), NF(Inf), x / (y + eps(NF)))

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

@kwdef struct QuadraticFunction{NF}
    "Quadratic coefficient"
    a::NF = 1.0
    "Linear coefficient"
    b::NF = 1.0
    "Constant coefficient"
    c::NF = 0.0
end

function (func::QuadraticFunction{NF})(x::NF) where {NF}
    (; a, b, c) = func
    y = a * x^2 + b * x + c
    return y
end
