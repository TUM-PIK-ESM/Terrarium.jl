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
deduplicate(values::Tuple) = Tuple(unique(values))
deduplicate(f, values::Tuple) = Tuple(unique(f, values))

"""
    $TYPEDSIGNATURES

Filters out all entries from `nt` that exist in `other`; like `setdiff` but for `NamedTuple`.
"""
@inline function ntdiff(nt::NamedTuple, other::NamedTuple{excluded}) where {excluded}
    nt_keys = keys(nt)
    filtered_keys = filter(∉(excluded), nt_keys)
    return NamedTuple{filtered_keys}(map(k -> getproperty(nt, k), filtered_keys))
end

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
