"""
    findfirst_z(i, j, condition_func, z_nodes, field)

2D kernel function that finds the first coordinate in `z_nodes` where `condition_func(field[i, j, k])`.
This implementation performs a linear scan over the z-axis and thus has time complexity O(N_z).
"""
@inline function findfirst_z(i, j, condition_func, z_nodes, field)
    i, j = idx
    n = length(z_nodes)
    idx = -1
    for k in 1:n
        if idx < 0 && condition_func(field[i, j, k])
            idx = k
        end
    end
    return idx > 0 ? z_nodes[idx] : z_nodes[n]
end
