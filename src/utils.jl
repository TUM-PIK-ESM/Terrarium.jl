"""
Alias for `Type{Val{x}}`
"""
const TypeVal{x} = Type{Val{x}} where {x}

"""
    $SIGNATURES

Concatenate one or more tuples together.
"""
tuplejoin() = tuple()
tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

"""
    $SIGNATURES
    
Filter out duplicates from the given tuple. Note that this method is not type stable or allocation-free!
"""
merge_duplicates(values::Tuple) = Tuple(unique(values))

"""
    merge_recursive(nt1::NamedTuple, nt2::NamedTuple)

Recursively merge two nested named tuples. This implementation is loosely based on the one in
[NamedTupleTools](https://github.com/JeffreySarnoff/NamedTupleTools.jl) authored by Jeffrey Sarnoff.
"""
merge_recursive(nt1::NamedTuple, nt2::NamedTuple, nts...) = merge_recursive(merge_recursive(nt1, nt2), nts...)
function merge_recursive(nt1::NamedTuple, nt2::NamedTuple)
    keys_union = keys(nt1) ∪ keys(nt2)
    pairs = map(Tuple(keys_union)) do key
        v1 = get(nt1, key, nothing)
        v2 = get(nt2, key, nothing)
        key => merge_recursive(v1, v2)
    end
    return (; pairs...)
end
# Recursive merge cases
merge_recursive(val, ::Nothing) = val
merge_recursive(::Nothing, val) = val
merge_recursive(_, val) = val # select second value to match behavior of merge

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
safediv(x::NF, y::NF) where {NF} = ifelse(iszero(y), Inf, x / (y + eps(NF)))

"""
    $SIGNATURES

Return a function `f(z)` that linearly interpolates between the given `knots`.
"""
function piecewise_linear(knots::Pair{<:LengthQuantity}...; extrapolation=Interpolations.Flat())
    # extract coordinates and strip units
    zs = collect(map(ustrip ∘ first, knots))
    ys = collect(map(last, knots))
    @assert issorted(zs, rev=true) "depths must be sorted in descending order"
    interp = Interpolations.interpolate((reverse(zs),), reverse(ys), Interpolations.Gridded(Interpolations.Linear()))
    return Interpolations.extrapolate(interp, extrapolation)
end

# fastmap and fastiterate

# Note that fastmap and fastiterate are borrowed (with self permission!) from CryoGrid.jl:
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

"""
    fastiterate(f!::F, iters::NTuple{N,Any}...) where {F,N}

Same as `fastmap` but simply invokes `f!` on each argument set without constructing a tuple.
"""
@generated function fastiterate(f!::F, iters::NTuple{N,Any}...) where {F,N}
    expr = Expr(:block)
    for j in 1:N
        push!(expr.args, :(f!($(map(i -> :(iters[$i][$j]), 1:length(iters))...))))
    end
    push!(expr.args, :(return nothing))
    return expr
end
fastiterate(f!::F, iters::NamedTuple) where {F} = fastiterate(f!, values(iters))
