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
