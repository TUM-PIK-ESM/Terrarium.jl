struct NumberFormatAdaptor{NF} end

"""
    Adapt.adapt(::NumberFormatAdaptor{NF}, obj) where {NF<:Number}

Adaptor that reconstructs arbitrary data structures with all numeric values
converted to the specified number format `NF`.
"""
function Adapt.adapt(::NumberFormatAdaptor{NF}, obj) where {NF<:Number}
    vals = map(NF, flatten(obj, flattenable, Number))
    return reconstruct(obj, vals, flattenable, Number)
end
