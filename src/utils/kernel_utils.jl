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
    min_zᵃᵃᶠ(i, j, k, grid, x)
    min_zᵃᵃᶠ(i, j, k, grid, f, args...)

Computes the field or function at the vertical (z-axis) face by taking the `min` of the two adjacent vertical layers.
"""
@inline min_zᵃᵃᶠ(i, j, k, grid, c) = @inbounds min(c[i, j, k], c[i, j, k - 1])
@inline min_zᵃᵃᶠ(i, j, k, grid, f, args...) = @inbounds min(f(i, j, k, grid, args...), f(i, j, k - 1, grid, args...))
