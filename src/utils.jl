"""
    $SIGNATURES

Concatenates one or more tuples together.
"""
tuplejoin() = tuple()
tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

"""
    $SIGNATURES
    
Filters out duplicates from the given tuple. Note that this method is not type stable or allocation-free!
"""
merge_duplicates(values::Tuple) = Tuple(unique(values))

"""
    $SIGNATURES

Evaluates `x / y` unless `iszero(y)` is true, then returns zero.
"""
safediv(x, y) = ifelse(iszero(y), zero(x), x / y)

"""
    $SIGNATURES

Returns a function `f(z)` that linearly interpolates between the given `knots`.
"""
function piecewise_linear(knots::Pair{<:LengthQuantity}...; extrapolation=Interpolations.Flat())
    # extract coordinates and strip units
    zs = collect(map(ustrip ∘ first, knots))
    ys = collect(map(last, knots))
    @assert issorted(zs, rev=true) "depths must be sorted in descending order"
    interp = Interpolations.interpolate((reverse(zs),), reverse(ys), Interpolations.Gridded(Interpolations.Linear()))
    return Interpolations.extrapolate(interp, extrapolation)
end

# fastmap

# Note that fastmap is borrowed (with self permission!) from CryoGrid.jl:
# https://github.com/CryoGrid/CryoGrid.jl/blob/master/src/Utils/Utils.jl
"""
    fastmap(f::F, iter::NTuple{N,Any}...) where {F,N}

Same as `map` for `NTuple`s but with guaranteed type stability. `fastmap` is a `@generated`
function which unrolls calls to `f` into a loop-free tuple construction expression.
"""
@generated function fastmap(f::F, iters::NTuple{N,Any}...) where {F,N}
    expr = Expr(:tuple)
    for j in 1:N
        push!(expr.args, :(f($(map(i -> :(iters[$i][$j]), 1:length(iters))...))))
    end
    return expr
end

"""
    fastmap(f::F, iter::NamedTuple...) where {F}

Same as `map` for `NamedTuple`s but with guaranteed type stability. `fastmap` is a `@generated`
function which unrolls calls to `f` into a loop-free tuple construction expression. All named
tuples must have the same keys but in no particular order. The returned `NamedTuple` 
"""
@generated function fastmap(f::F, nts::NamedTuple...) where {F}
    expr = Expr(:tuple)
    # get keys from first named tuple
    keys = nts[1].parameters[1]
    for key in keys
        push!(
            expr.args,
            :($key = f($(map(i -> :(nts[$i].$key), 1:length(nts))...)))
        )
    end
    return expr
end
