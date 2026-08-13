"""
    @assert_kernel cond [text]

Drop-in replacement for `Base.@assert` that is elided inside a Reactant compile context.
"""
macro assert_kernel(cond, text...)
    # Splice the `@assert` in with the caller's line number so failures point at the callsite
    # rather than at this macro definition.
    assertion = Expr(:macrocall, GlobalRef(Base, Symbol("@assert")), __source__, esc(cond), map(esc, text)...)
    return quote
        if !$(ReactantCore.within_compile)()
            $assertion
        end
        nothing
    end
end

struct KernelFunction{autonomous, Func, Args}
    func::Func
    args::Args
end

function (op::KernelFunction{false})(var, grid, clock, fields, args...)
    loc = map(typeof, location(vardims(var)))
    return KernelFunctionOperation{loc...}(op.func, get_field_grid(grid), clock, fields, args..., op.args...)
end

function (op::KernelFunction{true})(var, grid, clock, fields, args...)
    loc = map(typeof, location(vardims(var)))
    return KernelFunctionOperation{loc...}(op.func, get_field_grid(grid), fields, args..., op.args...)
end

"""
    kernel(func, args...; clock = false)

Return a `KernelFunction` that lazily constructs a `KernelFunctionOperation` from the given
`func` and tuple of `args` when invoked, i.e:

```juila
ctor = kerenl(my_function, arg1, arg2)
...
kfo = ctor(grid, clock, fields)
```

This is intended to be used as a constructor for `AuxiliaryVariable`s:

```julia
myvar(i, j, k, grid, fields) = clamp(fields.x[i, j, k], zero(eltype(grid)), one(eltype(grid)))

auxvar = auxiliary(:myvar, XYZ(), kernel(myvar))
```
"""
function kernel(func, args...; clock = false)
    return KernelFunction{!clock, typeof(func), typeof(args)}(func, args)
end

"""
    findfirst_z(i, j, condition_func, z_nodes, field)

2D kernel function that finds the first coordinate in `z_nodes` where `condition_func(field[i, j, k])`.
This implementation performs a linear scan over the z-axis and thus has time complexity O(N_z).
"""
@propagate_inbounds function findfirst_z(i, j, condition_func, z_nodes, field)
    n = length(z_nodes)
    idx = -1
    for k in 1:n
        found = (idx < 0) & condition_func(field[i, j, k])
        idx = ifelse(found, k, idx)
    end
    # Select the index (not the nodes): ifelse evaluates both branches, so indexing in the
    # unselected branch with idx = -1 would be out of bounds.
    return z_nodes[ifelse(idx > 0, idx, n)]
end

"""
    pad_indices(field, indices)

Expand the index tuple `indices` to the dimensionality of `field`, padding with trailing ones.

`Field`s are always three-dimensional and only define `setindex!` for exactly three indices, so a 2D
write falls back to `Base`'s generic `AbstractArray` path. That path calls `axes(::AbstractField)` →
`size(::AbstractGrid, ...)`, which is a dynamic dispatch on the device and fails to compile under
Reactant. Solvers operating on a 2D (`XY`) slice therefore pad their indices before writing to the
target field.
"""
@inline pad_indices(field, indices::NTuple{N, Integer}) where {N} = ntuple(d -> d <= N ? indices[d] : 1, Val(ndims(field)))

"""
    min_zᵃᵃᶠ(i, j, k, grid, x)
    min_zᵃᵃᶠ(i, j, k, grid, f, args...)

Computes the field or function at the vertical (z-axis) face by taking the `min` of the two adjacent vertical layers.
"""
@inline min_zᵃᵃᶠ(i, j, k, grid, c) = @inbounds min(c[i, j, k], c[i, j, k - 1])
@inline min_zᵃᵃᶠ(i, j, k, grid, f, args...) = @inbounds min(f(i, j, k, grid, args...), f(i, j, k - 1, grid, args...))
