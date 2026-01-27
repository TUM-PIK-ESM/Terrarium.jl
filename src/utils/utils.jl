"""
Alias for `Type{Val{x}}`
"""
const ValType{x} = Type{Val{x}} where {x}

"""
Alias for `Union{Nothing, T}` indicating that an argument or field of type `T` is optional and
can be replaced with `nothing`.
"""
const Optional{T} = Union{Nothing, T}

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

include("tuple_utils.jl")
include("interpolation_utils.jl")
include("kernel_utils.jl")
include("adaptors.jl")
